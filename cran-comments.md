## v0.1.2
* Now code is heavily documented and manual has been made more precise
* Made lik computation more efficient
* Made interpolation near origin default variogram fitting choice (linear is second, exp is third)
* Added option for user to specify xlims and ylims of BAUs 
* Removed timeline object (was redundant)
* Internally forced all time objects to be POSIXct

## v0.1.1 Since first attempt to upload on CRAN
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

