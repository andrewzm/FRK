Fixed Rank Kriging
================

This package is operational but still in its early stages of development. It is based on the work by Cressie and Johannesson (2008) and carries out using Fixed Rank Kriging (FRK) using large datasets. A special feature of this package is that it allows for a Hadoop backend, allowing for the construction of massive L3 products.

Package: FRK
Type: Package
Title: Fixed Rank Kriging
Version: 0.1.0
Date: 2015-05-15
Author: Andrew Zammit-Mangion
Maintainer: Andrew Zammit-Mangion <andrewzm@gmail.com>
VignetteBuilder: knitr
Description: Fixed Rank Kriging
BugReports: http://github.com/andrewzm/FRK/issues
Depends:
    R (>= 3.1)
Suggests:
    testthat,
    reshape2,
    shapefiles,
    knitr,
    Rhipe
Imports:
    methods,
    Matrix,
    digest,
    gstat,
    ggplot2,
    dplyr,
    scales,
    INLA,
    grid,
    fields,
    plyr,
    PBSmapping,
    akima,
    parallel,
    sp
License: GPL (>= 2)
NeedsCompilation: yes
LazyData: true
