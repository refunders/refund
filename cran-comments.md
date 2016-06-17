## Test environments
* local windows 8 x64, R 3.1.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or NOTEs.

There is 1 WARNING: "'library' or 'require' call not declared from: ‘dtw’" this is because the package is called via the "method="dtw", window.type="sakoechiba"" options to dist().
