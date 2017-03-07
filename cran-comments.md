## Since first attempt to upload on CRAN
* Native routines are now registered in FRK-init.c
* vignettes now do not include packages that are not on CRAN and re-building of vignettes should be successful
* all copyright is now acknowledged (checked using grep in all folders)
* examples now execute without randomly stopping

## Test environments
* local ubuntu 14.04, R 3.2.0
* local ubuntu 14.04, R 3.3.2
* local Win64, R 3.3.2
* win-builder, R 3.3.3
* win-builder, R-devel

## R CMD check results
There were no ERRORs, WARNINGs.

There were two NOTEs under win-builder, R 3.3.3
  - this is a new submission. 
  - the package installed size is 8.3Mb.

There is one additional NOTE under win-builder R-devel
  - this package enhances dggrids which is not available for checking.
  
## Downstream dependencies
There are no downstream dependencies

