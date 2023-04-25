## Resubmission
This is a resubmission. In this version I have:

* Add quotes around the title of the reference

* Add \value to .Rd files, specifically in file plot_estimates.Rd, where it was missing previously.

* Replace \dontrun{} with \donttest{} as those examples are not executable in < 5 sec.

* Replace cat() with message().


## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Downstream dependencies
There are currently no downstream dependencies for this package.
