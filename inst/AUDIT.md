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

## Résultats de l'étape 2 (A1, et une partie de A3)

**Activation en bloc** (`Optimizer::select_violators`, `working_set` dans `OptimizerSparse.h` et
`OptimizerGroup.h`) :

- à chaque itération externe, on active jusqu'à `maxadd` variables (ou groupes) inactives parmi
  les plus fortes violations des KKT, au lieu de la seule plus forte violation globale ;
- nouveau paramètre de contrôle `maxadd` : 10 par défaut pour les modèles parcimonieux et Lava,
  5 pour les modèles à groupes (valeurs choisies sur les mesures ci-dessous) ;
- l'activation est plafonnée pour ne pas dépasser `maxfeat + 1` variables actives ;
- le statut « max # of iterate reached » ne dépend plus du compteur d'itérations mais du gap
  final (un λ convergé exactement à la dernière itération n'est plus signalé en échec) ;
- correction au passage de l'initialisation QUADRA des groupes, qui lisait `grad(grp_in)` avec un
  numéro de groupe au lieu du gradient des variables du groupe.

Les *strong rules* n'ont pas été implémentées : dans une méthode d'ensemble actif, elles ne
restreignent pas les variables à activer (toutes les violatrices sont dans l'ensemble fort) et
servent seulement à éviter le test KKT sur toutes les variables, qui est en O(p) et n'est pas le
goulot. Les *Gap Safe rules* restent une piste.

**Solveurs proximaux** (`Optimizer.cpp`). L'activation en bloc a révélé une fragilité existante :
sur des données mal conditionnées (prostate non normalisée, normes de colonnes de 4,6 à 633),
PGD et FISTA s'arrêtaient après une seule itération interne, car le critère ‖x⁺ − x‖ < 1e-7 est
proportionnel au pas 1/L. La référence elle-même ne convergeait pas sur ce cas et ne passait le
test (`test-enet-reference.R`, tolérance 1e-2) que par chance. Deux corrections :

- critère d'arrêt sur la norme de l'application gradient, L·‖x⁺ − x‖, invariante au pas ;
- garde-fou pour Anderson (PGD) : si le résidu augmente après un pas extrapolé, ce pas est rejeté,
  on revient au pas proximal simple et l'historique est vidé.

Sur 50 tirages de `lambda2` (le test tire `lambda2` au hasard sans graine), PGD avec activation en
bloc échouait 38 à 50 fois avant ces corrections et 0 fois après ; sur la graine 1, pgd et fista
atteignent l'objectif de référence à 3e-14 près (1e-5 dans la version d'origine).

**Validation** : 128/128 tests. Là où les deux versions convergent, coefficients identiques à la
référence à 1e-14 (lasso, elastic-net), 1e-6 (MCP, tolérance de la LLA), 1e-8 (FISTA). Pour le
group-lasso fista/pgd, deux groupes supplémentaires à coefficients ≤ 1,8e-4 au dernier λ, avec un
objectif légèrement *plus bas* que la référence (736,782724 contre 736,782733).

**Choix de `maxadd`** (secondes, même protocole) :

| Cas | 1 | 5 | 10 | 20 | 100 |
|---|---|---|---|---|---|
| Elastic-net λ₂=1, n=500, p=10 000, quadra | 34,1 | 12,8 | 9,4 | 7,9 | 34,3 |
| Elastic-net λ₂=1, n=200, p=2 000, quadra | 1,06 | 0,65 | 0,58 | 0,61 | 0,97 |
| Lasso n=500, p=10 000, quadra | 1,11 | 0,79 | — | 0,80 | 0,93 |
| Group-lasso pgd | 0,91 | 0,56 | — | 1,24 | 1,35 |

Au-delà d'une vingtaine de variables, le solveur de Newton passe son temps à retirer les variables
ajoutées en trop.

**Temps avec les valeurs par défaut** (secondes, `b1dfd64` → branche) :

| Cas | quadra | fista | pgd |
|---|---|---|---|
| Elastic-net λ₂=1, n=500, p=10 000 | 112,4 → **7,2** (non convergés : 28 → 1) | — | — |
| Lasso n=500, p=10 000 | 2,05 → **0,77** | 1,65 → 0,94 | 1,18 → 0,78 |
| Lasso n=200, p=2 000 | 0,08 → 0,08 | 0,14 → 0,10 | 0,13 → 0,11 |
| Lasso n=2 000, p=5 000 (k ≤ 41) | 0,87 → 0,97 | 0,98 → 1,09 | 0,85 → 0,95 |
| Group-lasso n=300, p=3 000 | 2,18 → 2,60 (toujours 25 non convergés, cf. A2) | 2,77 → 2,27 | 1,45 → 0,97 |
| Lava p=2 000 | 0,43 → 0,42 | — | — |

Le lasso p=10 000 est désormais à 2× de glmnet (0,39 s) au lieu de 5×. Seul point négatif : un
ralentissement d'environ 10 % sur le lasso n=2 000 avec un très petit ensemble actif (k ≤ 41).

Sur l'elastic-net, le chemin s'arrête à 84 λ au lieu de 89 : la version convergée atteint
réellement `maxfeat` plus tôt, alors que la référence, non convergée, y arrivait en retard.

## Résultats de l'étape 3 (A3)

Les deux premiers points de A3 (critère d'arrêt invariant au pas, garde-fou Anderson) ont été
faits à l'étape 2. Cette étape termine le sujet des solveurs proximaux (`Optimizer.cpp`,
`working_set`).

- **Redémarrage adaptatif de FISTA** (O'Donoghue & Candès, 2015) : le moment est remis à zéro
  quand il pointe contre l'application gradient, soit (y_k − x_{k+1})ᵀ(x_{k+1} − x_k) > 0.
- **Estimation de L par Lanczos** au lieu de la puissance itérée. Mesures sur des matrices de Gram
  d'ensembles actifs croissants (n=500, k jusqu'à 2 000) :
  - données corrélées (ρ=0,3) : la puissance itérée converge en ~5 itérations ; Lanczos est exact
    à 1e-16 en 10 pas ;
  - données non corrélées (ρ=0), valeurs propres du haut serrées : la puissance itérée ne converge
    pas en 50 itérations et **sous-estime λmax jusqu'à 2,9 %** (8 % avec les anciens réglages,
    15 itérations), au-delà de la marge de 1 %, d'où un pas trop grand ; la borne de Gershgorin est
    2,9 fois trop grande. Avec 30 pas de Lanczos, la borne θ + résidu est toujours ≥ λmax et au
    plus 0,5 % au-dessus.

  Implémentation : au plus 30 pas, arrêt dès que le résidu de Ritz est ≤ 1e-3 θ, L = max(θ +
  résidu, 1,01 θ, max diag). Départ déterministe : l'ancienne initialisation par `randu`
  consommait le générateur aléatoire de R et modifiait la suite aléatoire de l'utilisateur.
- **Réestimation de L après un retrait** : `set_changed` était remis à `false` au début d'une
  itération sans ajout, si bien qu'un retrait de variable à l'itération précédente ne déclenchait
  pas de réestimation (L périmé mais majorant, donc pas dangereux, seulement trop prudent).

**Validation** : 128/128 tests ; 0 échec sur les 50 tirages prostate en pgd et fista ; supports
identiques à l'étape 2 et coefficients à ≤ 1,4e-7 près (fista), ≤ 5e-9 (pgd).

**Temps** (secondes, étape 2 → étape 3 ; itérations internes entre parenthèses) :

| Cas | fista | pgd |
|---|---|---|
| Group-lasso n=300, p=3 000, ρ=0,3 | 2,20 → **0,62** (40 569 → 7 303) | 1,01 → 1,00 |
| Group-lasso, ρ=0 | 0,65 → 0,44 (7 371 → 2 857) | 0,49 → 0,52 |
| Lasso n=500, p=10 000, ρ=0,3 | 1,02 → 0,85 (45 025 → 9 193) | 0,87 → 0,87 |
| Lasso n=500, p=10 000, ρ=0 | 1,04 → 0,87 (21 104 → 6 054) | 0,97 → 0,98 |
| Lasso n=200, p=2 000 | 0,10 → 0,06 (29 835 → 6 753) | 0,11 → 0,10 |
| `bounded_reg` n=300, p=1 000, λ₂=1 | 11,06 → **2,76** (64 238 → 13 556) | — |

Sur le lasso p=10 000, FISTA fait 5 fois moins d'itérations internes mais ne gagne que 17 % :
le coût fixe par itération externe domine désormais (voir P5).

**Non fait** : tolérance interne adaptée au gap externe (résolution inexacte). Gain attendu
faible maintenant que le coût des itérations internes n'est plus dominant.

## Résultats de l'étape 4 (A2)

**Nouveau solveur QUADRA pour les modèles à groupes** (`GroupOptimizer::quadratic` et
`block_solve` dans `OptimizerGroup.h`) : descente par blocs exacte sur les groupes actifs.

- **Test de nullité exact** pour chaque groupe : β_g = 0 si la norme duale du résidu partiel
  (soft-thresholdé si α > 0) est ≤ λw. Il vaut pour les trois normes (`penalty_.optimality`).
- **Group-lasso (α = 0)** : solution exacte du bloc. Avec l'EVD H_g = V diag(d) Vᵀ déjà en
  cache, t = ‖β_g‖ résout Σ c_i² / (d_i t + μ)² = 1 avec c = Vᵀr ; la fonction est convexe et
  décroissante, donc Newton parti de la gauche converge de façon monotone.
- **Pas de Newton global** après chaque passe, sur les groupes non nuls, avec le hessien exact
  H + Σ_g μ_g/t_g (I − u_g u_gᵀ) et une recherche linéaire d'Armijo. Sans lui, la descente par
  blocs plafonnait à 1 000 passes sur les données du test (n = 50, p = 95, 5 groupes corrélés).
- **Sparse-group, coop, l1/l∞** : FISTA avec redémarrage sur le bloc, avec la prox exacte de la
  pénalité et L = max de l'EVD du groupe.
- Les groupes nuls sont retirés à la fin, une fois la descente convergée ; les nouveaux groupes
  partent de 0 ; l'EVD est maintenue dès que la méthode est quadra (`factmat` dans
  `group_sparse_lm` et `group_lava`). Le poids du terme L1 est désormais celui de la pénalité
  (non pondéré, comme dans `elt_norm` et `proximal`).

**Bug corrigé dans la prox l1/l∞** (`PenaltyGroup.cpp`) : quand ‖x_g‖₁ ≤ λw, la prox de
λw‖·‖∞ vaut 0, mais le code laissait x_g inchangé. FISTA et PGD en l1/l∞ convergeaient vers une
mauvaise solution.

**Validation** : 128/128 tests, dont le test de temps « quadra plus rapide que fista »
(0,04 s contre 0,10 s). Objectifs pénalisés calculés indépendamment en R : quadra, fista et pgd
coïncident à 4e-12 près sur tous les λ, pour les quatre pénalités.

**Temps** (n=300, p=3 000, groupes de 10, λ₂ = 0 ; secondes) :

| Pénalité | quadra avant | quadra après | fista après | grpreg |
|---|---|---|---|---|
| Group-lasso, ρ=0,3 | 1,73 (22 λ non convergés) | **0,37** | 0,40 | 0,18 |
| Group-lasso, ρ=0 | 1,66 (21 non convergés) | 0,48 | 0,40 | 0,14 |
| Sparse-group α=0,5 | 24,5 (99 non convergés, objectif ×25) | **0,43** | 0,31 | — |
| Coop | 1,40 (19 non convergés) | 0,74 | 0,43 | — |
| l1/l∞ | inutilisable | 0,72 | 57,0 → **0,58** (bug de prox) | — |

« Non convergés » compte les λ au statut autre que `converged`. Après correction, le seul λ
restant dans chaque cas est l'arrêt du chemin sur « max # of feature reached » (570 variables
actives pour `maxfeat` = 600), identique en fista : tous les λ convergent.

## Résultats de l'étape 5 (P3, P5, P6)

**P3 — Downdate de Givens** (`ActiveSet.h`) : rotations appliquées sur place sur les lignes
(k, k+1), sans temporaires Armadillo. Micro-benchmark (suppression d'une colonne au quart) :
0,13 → 0,04 ms pour k = 200, 7,2 → 5,4 ms pour k = 1 000, 32 → 26 ms pour k = 2 000 (les deux
`shed_*` restent). Résultat identique à 1e-13.

**P5 — Boucles proximales** (`working_set` des deux optimiseurs) : poids, tailles de groupes et
X'y restreints à l'ensemble actif calculés une fois par appel du solveur au lieu d'une fois par
itération interne. Gain faible (group-lasso fista 0,39 → 0,34 s) : le passage de `std::function`
à un paramètre template n'a pas été fait, le coût dominant étant ailleurs (voir ci-dessous).

**P6 — Points secondaires faits** :

- `criteria()` (R) évaluait l'active binding `deviance` 5 fois, chacune refaisant
  `X %*% coef`, `sweep` et `apply` : il est évalué une fois ; `residuals` est vectorisé. Lasso
  n=500, p=10 000, temps total R : 0,73 → 0,58 s.
- Ridge et Lava avec une structure S diagonale à diagonale positive (le cas par défaut) :
  C⁻¹ = diag(1/√s) appliqué comme mise à l'échelle des colonnes, sans Cholesky ni inverse
  p × p ; côté R, `CholStruct()` n'est appelé que si S n'est pas de ce type. n=200 :

  | p | ridge avant → après | lava avant → après |
  |---|---|---|
  | 4 000 | 0,61 → 0,21 s | 1,68 → 0,33 s |
  | 8 000 | 1,98 → 0,40 s | 11,6 → **0,63 s** |

  Coefficients identiques à 4e-15, degrés de liberté à 9e-14.
- `BoundedRegression` : coefficients pré-alloués (plus de `join_rows` dans la boucle) ; `get_df`
  en `accu(SUU % C)` avec une seule passe sur les non-zéros de S ; matrice de Gram en forme A'A
  (`syrk`). `bounded_reg` n=300, p=1 000 : 1,18 → 1,07 s.
- `RidgeRegression` : coefficients pré-alloués, produit sans `diagmat`.
- `wrapper_FusedLasso.cpp` : normalisation des colonnes par itérateur sur les non-zéros.
- `src/Makevars` et `src/Makevars.win` : les objets dépendent des en-têtes ;
  `-DARMA_NO_DEBUG` ajouté sous Windows.
- Validation croisée : les données de chaque pli sont construites dans la tâche du pli
  (`DataModel$splitFold`) au lieu de construire les K plis à l'avance ; même temps, K copies de
  X en moins en mémoire simultanément. Erreurs de CV identiques à 2e-15.

**Validation** : 128/128 tests ; résultats identiques à l'étape 4 (lasso, group-lasso, ridge et
lava à ≤ 4e-15, `bounded_reg` à 3,7e-7, critères d'information à 9e-13).

**P6 non fait** : conversion de X dense en creux et normalisation non centrée ni pondérée dans
FusedLasso (changement de comportement, à décider) ; warm start entre valeurs de λ₂ et passage
d'indices au C++ pour la validation croisée.

**Où part le temps maintenant.** Lasso n=500, p=10 000 : environ 0,4 s en C++ quel que soit le
solveur, dont 90 % sur les λ où des variables entrent. Le coût est le calcul des nouvelles
colonnes X'WX_j sur les p lignes (O(np) par variable), imposé par le stockage de `XTXA_` : c'est
l'objet de P4.

## Feuille de route

| Étape | Contenu | Risque | Statut |
|---|---|---|---|
| 1 | P1 + P2 (ensemble actif pré-alloué, solves triangulaires) | faible, couvert par les tests | fait (33f7102) |
| 2 | A1 (activation en bloc) + statut de convergence ; critère d'arrêt et garde-fou Anderson (A3) | moyen | fait (537c847, 3a28cd8) |
| 3 | A2 (descente par blocs exacte + Newton pour les groupes) ; bug de prox l1/l∞ | moyen | fait |
| 4 | A3 restant (restart FISTA, estimation de L par Lanczos) | faible | fait (934346e) |
| 5 | P3, P5, P6 | faible | fait |
| 5 bis | P4 (gradient par le résidu, mémoire O(pk)) | moyen | à faire |
| 6 | A4 (CD + working set) | élevé | à évaluer |
