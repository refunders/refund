## CRAN submission 0.1-30

Version: 0.1-29
Check: S3 generic/method consistency

## CRAN submission 0.1-29

* Fixed errors due to failure to find CanadianWeather data in examples

## CRAN submission 0.1-28

The Date field is not in ISO 8601 yyyy-mm-dd format.

* Fixed for resubmission



## Test environments
* local windows 8 x64, R 3.1.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or NOTEs.

There is 1 WARNING: "'library' or 'require' call not declared from: ‘dtw’" this is because the package is called via the "method="dtw", window.type="sakoechiba"" options to dist().
