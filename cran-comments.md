## v0.1.7
+ Made compatible with new update of Hmisc
+ Added support for change of support spatio-temporal prediction
+ Changed from pred_locs to newdata in arguments to SRE.predict

## v0.1.6
* Added scale_aperture to have control on aperture of basis functions
* Adde remove_basis to easily remove basis functions
* Made FRK compatible with Rcpp 0.12.12

## v0.1.4
* Robust basis function placement in 2D
* Added offset parameter when constucting mini BAUs to avoid
  underflow
* Fixed problem of estimating length scale from one data point
* Clarified copyright and license information in files

## v0.1.3
* Important added functionality for using package with only point data and point predictions, added on request.

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

