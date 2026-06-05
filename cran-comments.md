## quadrupen 1.0-0	(2026-06-05)

- major updates
  - complete rewriting of R code using R6 classes
  - complete rewriting of C++ code using template and OO style programming
  - included 'FusedLasso' from the archived package by Holger Hoefling (fixed CRAN's complaints)
  - added group-lasso/group-elastic and variants (group penalty: l1/l2, l1/linf, cooperative Lasso)
  - added sparse group-lasso/group-elastic and variant (group penalty: l1/l2, l1/linf, cooperative Lasso)
  - added lava and/post-lava (combination of sparse and dense regularization, Chernozukov et al, 2017)
  - extended to group-lava (group penalty: l1/l2, l1/linf, cooperative Lasso)
  - added mcp and scad (concave penalties, Chernozukov et al, 2017)
  - added refit version of Lasso/Elastic-Net ("relaxed" Lasso/Enet)
  - changing many parameters (badly) named, do not expect backward compatibility
  - added vignettes
- minor updates
  - Integration of changes from CRAN versions from 0.2-4 to 0.2-13
  - set up github workflow for pkgdown page
  - various fixes, more testing

## Tested environments

* tested locally on Ubuntu Linux 24.04.4 LTS, R-release, GCC

* tested remotely with github-action

- Linux ubuntu 24.04, R-release
- Linux ubuntu 24.04, R-oldrel
- Linux ubuntu 24.04, R-devel
- Windows Server 2025, R-release, 64 bit
- macOS 15, R-release

* tested remotely with win-builder (R version 4.5.3, R unstable, R version 4.6.0)

all status OK except 1 NOTE (False Positive, MCP is the name of a penalty, not a misspelled word)

Possibly misspelled words in DESCRIPTION:
  MCP (9:165)
  

## Local R CMD check results

── R CMD check results ────────────────────────────────────────────────────────────────────────────────── quadrupen 1.0-0 ────
Duration: 4m 23.1s

0 errors | 0 warnings | 0 note
