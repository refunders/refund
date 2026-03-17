## CRAN submission 0.1-40

### Breaking changes

* Default `sandwich` option for `pffr()` changed from `"none"` to `"cluster"`.
  An informational `message()` is printed when `sandwich` is not explicitly
  supplied.
* `pffrGLS()` / `pffr_gls()` now error with a deprecation message. Users are
  directed to `pffr()` with `sandwich = "cluster"` or `sandwich = "cl2"`.

### Major changes

* `pffr()` internal refactor into prepare-fit-postprocess pipeline.
* Function renames with deprecated wrappers: `pffrSim` -> `pffr_simulate`,
  `coefboot.pffr` -> `pffr_coefboot`, `qq.pffr` -> `pffr_qq`,
  `pffr.check` -> `pffr_check`.
* New CL2 leverage-adjusted cluster-robust covariance option.
* Improved AR(1) support with automatic `bam`/`fREML` switching.

## Test environments

* local Linux Mint 22.1 (R 4.5.2)
* R-hub v2 (Linux, Windows, macOS)

## R CMD check results

0 errors | 0 warnings | 1 note

* `checking for future file timestamps`: unable to verify current time
  (network issue, not a package problem)
