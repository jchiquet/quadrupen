# Audit de performance et d'algorithmes — quadrupen

- **Date** : 2026-09-16
- **Révision auditée** : `b1dfd64` (master, version 1.0-0)
- **Périmètre** : code C++ (`src/Quadrupen`, `src/FusedLasso`, wrappers) et orchestration R
  (validation croisée, sélection de stabilité), sous l'angle (i) performance de calcul pure,
  (ii) choix algorithmiques.

## Synthèse

Le code est propre et plusieurs optimisations sont déjà en place (cache de la constante de
Lipschitz, spécialisation `sp_mat` de la standardisation, `nth_element` dans FusedLasso).
Deux constats dominent :

1. **Le goulot principal n'est pas l'algèbre linéaire mais la gestion mémoire de l'ensemble
   actif** : `XTXA_` (p × k) est réalloué et recopié à chaque ajout ou retrait de variable, et
   les solves triangulaires matérialisent des copies k × k. Sur les chemins où l'ensemble actif
   devient grand, ces copies coûtent 10 à 40 fois le calcul utile.
2. **Le solveur QUADRA du group-lasso est approché** (point fixe tronqué, soft-threshold a
   posteriori) : il est à la fois lent et imprécis, et de nombreux λ ne convergent pas.

## Mesures

**Protocole.** Machine : Intel i7-12800H (6 cœurs performance avec hyperthreading, CPU 0–11 ;
8 cœurs efficacité, CPU 12–19). Chaque cas tourne dans un processus R séparé, BLAS/OpenMP à
1 thread (`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1`), épinglé par `taskset` sur un cœur
performance physique distinct (CPU 0, 2, 4, 6, 8), les deux versions d'un même cas s'exécutant
l'une après l'autre sur le même cœur. Sans épinglage, les tâches parallèles atterrissent sur des
cœurs efficacité ou sur des threads jumeaux et les temps sont faussés (jusqu'à ×2 sur les
témoins glmnet et grpreg). Données gaussiennes équi-corrélées (ρ = 0.3), 20 coefficients non
nuls, `set.seed(1)` par cas, 100 valeurs de λ. Temps en secondes, version `b1dfd64`.

| Cas | quadra | fista | pgd | Référence |
|---|---|---|---|---|
| Lasso n=200, p=2 000 | 0,07 | 0,12 | 0,13 | glmnet 0,02 |
| Lasso n=500, p=10 000 (k ≤ 249) | 2,25 | 2,00 | 1,42 | glmnet 0,55 |
| Lasso n=2 000, p=5 000 (k ≤ 41) | 0,88 | 0,99 | 0,89 | glmnet 0,76 |
| Elastic-net λ₂=1, n=500, p=10 000 (k ≤ 1 998) | 109,4, **28 λ non convergés** | — | — | — |
| Group-lasso n=300, p=3 000, groupes de 10 | 2,13, **25 non convergés, gap max 0,73** | 2,92 | 1,44 | grpreg 0,17 |
| Lava n=200, p=500 / 1 000 / 2 000 | 0,07 / 0,13 / 0,43 | — | — | — |

Lava p=300, n=250 / 500 / 1 000 : 0,15 / 0,18 / 0,49 s. Le calcul `trace(K * Proj_)` de
`get_df` n'est pas un problème : Armadillo évalue `trace(A*B)` sans former le produit.

Sur l'elastic-net, `maxiter = 5000` au lieu de 50 fait converger tous les λ sauf un, sans gain de
temps : il y a deux problèmes distincts, le coût par itération (P1, P2) et la limite
d'itérations externes (A1).

Micro-benchmark des opérations élémentaires de l'ensemble actif (n=500, p=10 000, en ms,
mesure non épinglée, ordre de grandeur) :

| k | réallocation `XTXA_` | X'Wx_j (calcul utile) | `XTXA_ * beta` | solve Cholesky (`trimatu(R).t()`) | 2 solves dans `quadratic` | `shed` R + XTXA |
|---|---|---|---|---|---|---|
| 500 | 13,6 | 1,3 | 1,5 | 0,4 | 0,7 | 28,6 |
| 1 000 | 27,1 | 1,7 | 2,6 | 5,2 | 3,1 | 59,7 |
| 2 000 | 54,0 | 1,3 | 5,4 | 25,0 | 16,2 | 134,7 |

Solve de Cholesky (avant + arrière), k = 200 / 1 000 / 2 000, en ms : `solve(trimat…)` actuel
0,12 / 5,9 / 21,1 ; substitutions écrites à la main 0,03 / 0,43 / 2,4 ; LAPACK `dtrtrs` 0,008 /
0,12 / 0,78.

## 1. Performance pure

### P1 — Réallocation de l'ensemble actif (gain attendu : ×5 à ×20 sur les grands ensembles actifs)

`src/Quadrupen/ActiveSet.h` (`add_var`, `add_vars`, `del_var`, `del_vars`)

- Chaque ajout réalloue `XTXA_` (p × k), `XATXA_` et `R_` puis recopie l'existant.
- Chaque retrait fait trois `shed_*` (recopies), et `del_vars` les répète variable par variable.

**Correctif** : stocker `XTXA_` dans un tampon à capacité croissante (croissance géométrique,
bornée par p) avec une taille logique k ; un retrait décale les colonnes en place ; un retrait
multiple compacte le tampon en une seule passe. Les produits `XTXA_ * v` se font sur une vue
sans copie.

### P2 — Solves triangulaires (gain attendu : ×5 sur ces opérations)

`src/Quadrupen/ActiveSet.h` (`update_Cholesky`, `solve_Gram`), `src/Quadrupen/OptimizerSparse.h` (`quadratic`)

`trimatu(R_).t()` / `trimatl(R_.t())` matérialisent une copie k × k à chaque appel ; dans
`quadratic` il manque en outre `solve_opts::fast` (estimation du conditionnement à chaque appel).

**Correctif** : appeler LAPACK `dtrtrs` (`arma::lapack::trtrs`) directement sur la mémoire de
`R_`, sans transposée ni copie : ×27 par rapport à l'existant et ×3 par rapport à des
substitutions écrites à la main (voir le micro-benchmark ci-dessus).

### P3 — Downdate de Givens

`src/Quadrupen/ActiveSet.h` (`downdate_Cholesky`)

Chaque rotation crée des temporaires Armadillo (`mat G`, `submat`, produit 2 × (p−k)).
Une boucle scalaire en place (ou `drot`) garde la complexité O(k²) avec une constante bien
plus faible.

### P4 — Mémoire O(pk) de `XTXA_`

Pour p = 10⁵ et k = 10³ : 800 Mo. Alternative : gradient par le résidu, g = X'(w ⊙ r), en
O(np) ou O(nnz) sans stockage. Stratégie hybride selon k·p vs n·p, et résidu d'office en creux.

### P5 — Allocations dans les boucles proximales

`src/Quadrupen/OptimizerSparse.h`, `src/Quadrupen/OptimizerGroup.h` (lambdas `prox`)

`weights.elem(set.A_)` et `grp_sizes_(G_)` sont recalculés à chaque itération interne, via
`std::function`, et la prox renvoie un `vec` neuf. Sortir ces quantités de la boucle, prox en
place, prox passée en paramètre template.

### P6 — Points secondaires

- `BoundedRegression.cpp`, `RidgeRegression.cpp` : `coef_ = join_rows(coef_, …)` dans la
  boucle sur λ, soit O(p·L²). Pré-allouer p × L.
- `BaseLava.h` (`lava_preprocess`), `RidgeRegression.cpp` : avec la structure par défaut
  S = Id, on calcule `chol(S)`, son inverse via `eye(p)` puis `Xwc * C_inv` : O(p³) + O(np²)
  et une matrice p × p dense inutiles. Détecter S diagonale et remplacer par une mise à
  l'échelle des colonnes. (Invisible jusqu'à p = 2 000, non mesuré au-delà.)
- `BoundedRegression::get_df` : `trace(SUU*C)` → `accu(SUU % C)` ; éviter la boucle
  `S_.at(i,j)` en O(k² log nnz).
- `RegressionData::precompute_XTX` : `X.t()*WX` passe par `gemm` ; former √w·X puis
  `X.t()*X` permet `syrk` (≈ ×2).
- `wrapper_FusedLasso.cpp` : `sp_x.col(i) /= normx(i)` probablement en O(p·nnz) (boucle par
  itérateur à la place) ; X dense convertie en creux ; normalisation ni centrée ni pondérée,
  contrairement au reste du package.
- `src/Makevars.win` : `-DARMA_NO_DEBUG` absent (vérifications de bornes actives sous Windows).
- `src/Makevars` : les objets ne dépendent pas des en-têtes. Presque tout le code est dans des
  `.h` templates, donc modifier un en-tête puis relancer `R CMD INSTALL` dans l'arborescence
  source réutilise des `.o` périmés. Ajouter une règle du type
  `$(OBJECTS): $(wildcard Quadrupen/*.h FusedLasso/*.h)`.
- Validation croisée (`R/QuadrupenFit-R6Class.R`, `cross_validate`) : X recopiée 2K fois en R
  puis encore une fois par fit en C++ ; pas de warm start entre valeurs de λ₂. Passer des
  indices ou des poids d'observation au C++ éviterait les copies.

## 2. Algorithmes

### A1 — Une seule variable activée par itération externe, `maxiter = 50`

`src/Quadrupen/OptimizerSparse.h` (`working_set`), `R/utils.R` (`optim_enet_default`)

Dès que plus de 50 variables entrent entre deux λ successifs, l'optimisation s'arrête en
« non convergé » (cause des 28 échecs de l'elastic-net mesuré). Pistes :

1. activer en bloc les m plus gros violateurs KKT (`add_vars` / `update_Cholesky_block`
   existent déjà) ;
2. pré-filtrer par la *sequential strong rule* (|g_j| > 2λ_{k+1} − λ_k) puis vérifier les KKT
   sur le reste ;
3. optionnellement, *Gap Safe screening rules* (Fercoq, Gramfort & Salmon, 2015), en s'appuyant
   sur le gap dual déjà calculé par `optimality_violation`.

### A2 — Group-lasso QUADRA approché

`src/Quadrupen/OptimizerGroup.h` (`quadratic`)

- une seule étape de point fixe μ = λw/‖β_g‖ par passe, au plus 15 passes, tolérance 1e-4 ;
- soft-threshold L1 appliqué après coup avec une courbure diagonale approchée ;
- test de nullité ‖β_g‖ < 1e-10 quasiment jamais atteint par ce point fixe.

**Remplacement** : test exact ‖r_g‖ ≤ λw_g ⇒ β_g = 0 ; sinon équation séculaire 1-D en
t = ‖β_g‖ résolue par quelques pas de Newton avec l'EVD déjà en cache (solution exacte du bloc).
Alternative à la grpreg, qui explique l'essentiel de l'écart ×200 : orthonormaliser chaque groupe
une fois (X_g = Q_g R_g), mise à jour de bloc en group soft-threshold fermé, descente par blocs
avec mise à jour du résidu en O(n·|g|).

### A3 — Solveurs proximaux

`src/Quadrupen/Optimizer.cpp`

- **FISTA sans restart** : redémarrage adaptatif (O'Donoghue & Candès, 2015), typiquement
  ×3 à ×10 (écart observé : 480 000 itérations FISTA contre 33 000 pour PGD + Anderson).
- **Warm start de Lipschitz inopérant** : réutilisé seulement à taille identique, alors qu'on ne
  réestime que lorsque l'ensemble change. Compléter l'ancien vecteur propre par 0 (ajout) ou
  retirer les entrées (retrait).
- **Estimation de L risquée** : 15 itérations de puissance + 1 % de marge peuvent sous-estimer L
  et faire diverger FISTA. Backtracking, ou borne de Gershgorin comme borne sûre.
- **Anderson (PGD) sans garde-fou** : n'accepter le pas extrapolé que s'il fait décroître
  l'objectif ou le résidu.
- **Tolérances internes absolues** (1e-7, 1e-9) : les rendre relatives et adaptées au gap externe.

### A4 — Direction « état de l'art » pour lasso / MCP / SCAD

Working set + solveur interne en descente par coordonnées accélérée par Anderson (celer :
Massias et al., 2018 ; skglm : Bertrand et al., 2022). Coût par itération en O(nnz), sans
factorisation. QUADRA resterait l'option haute précision pour de petits ensembles actifs.

## Anomalies relevées au passage (à vérifier)

- `src/Quadrupen/OptimizerGroup.h` (`working_set`, activation d'un groupe) : `grad(grp_in)`
  indexe le gradient **par variable** avec un numéro de **groupe** ; le signe
  d'initialisation est donc faux.
- `src/Quadrupen/BoundedRegression.cpp` (`solution_path`) : `sum(penalty_.optimality(...))`
  vaut ‖g‖₁,w − p·λ et non ‖g‖₁,w − λ ; le gap est très sous-estimé et la boucle s'arrête
  probablement dès la première itération.

## Résultats de l'étape 1 (P1 + P2)

**Modifications** (`src/Quadrupen/ActiveSet.h`, appelants dans `OptimizerSparse.h`,
`OptimizerGroup.h`, `SparseRegularizer.h`) :

- `XTXA_` remplacé par un tampon privé p × capacité à croissance géométrique (bornée par p) ;
  les produits passent par `XTXA_times(v)`, qui multiplie une vue sans copie ; un retrait
  compacte le tampon en place, en une seule passe pour `del_vars`.
- Solves triangulaires par LAPACK `dtrtrs` directement sur la mémoire de `R_`
  (`update_Cholesky`, `update_Cholesky_block`, `solve_Gram`, `inverse_Gram`) ; `quadratic`
  passe par `solve_Gram`.
- Si le facteur est dégénéré (pivot nul ou solution non finie), repli sur `solve(XATXA_, b)` ou
  `inv_sympd(..., allow_approx)`, comme le faisait le `solve` d'Armadillo sans `fast`.
- `XATXA_` et `R_` (k × k) restent réalloués à chaque ajout : leur coût est désormais du même
  ordre que le calcul utile.

**Validation** : 128/128 tests ; coefficients, intercepts et degrés de liberté identiques à
l'ancienne version à 1e-13 près sur lasso, MCP, elastic-net structuré avec refit, group-lasso
(quadra et pgd) ; nombres d'itérations internes et externes identiques ; résultats identiques sur
un design avec colonnes colinéaires.

**Temps** (même protocole, en secondes) :

| Cas | avant | après | gain |
|---|---|---|---|
| Elastic-net λ₂=1, n=500, p=10 000, quadra | 109,4 | 26,1 | ×4,2 |
| Lasso n=500, p=10 000, quadra | 2,25 | 1,01 | ×2,2 |
| Lasso n=500, p=10 000, pgd | 1,42 | 0,97 | ×1,5 |
| Lasso n=500, p=10 000, fista | 2,00 | 1,61 | ×1,2 |
| Lasso n=2 000, p=5 000 (k ≤ 41) | 0,88 | 0,92 | = |
| Lasso n=200, p=2 000 | 0,07 | 0,06 | = |
| Group-lasso quadra / fista / pgd | 2,13 / 2,92 / 1,44 | 1,98 / 2,43 / 1,03 | = / ×1,2 / ×1,4 |
| Lava p=2 000 | 0,43 | 0,41 | = |

Le gain croît avec la taille de l'ensemble actif, comme attendu. L'écart restant avec glmnet
(lasso p=10 000 : 1,01 s contre 0,42 s) et grpreg relève maintenant surtout des algorithmes
(A1, A2, A4).

## Feuille de route

| Étape | Contenu | Risque | Statut |
|---|---|---|---|
| 1 | P1 + P2 (ensemble actif pré-alloué, solves triangulaires) | faible, couvert par les tests | fait (non commité) |
| 2 | A1 (activation en bloc, strong rules) + révision de `maxiter` | moyen | à faire |
| 3 | A2 (bloc exact pour le group-lasso) | moyen | à faire |
| 4 | A3 (restart FISTA, warm start Lipschitz, garde-fou Anderson) | faible | à faire |
| 5 | P3–P6 | faible | à faire |
| 6 | A4 (CD + working set) | élevé | à évaluer |
