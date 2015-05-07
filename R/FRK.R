#' Multi-variate spatio-temporal modelling
#'
#' This package provides a set of tools for modelling in space and time. Currently focus is placed on the SPDE approach of Finn Lindgren et al. 2011
#' but extensions are straightforward to implement please contribute by forking this package and asking a pull request to include your work. The key
#' use of this package is to cater for multiple ST processes and observation processes with ease. Then a "compress" function can be called to decompose
#' the multi-variate ST fields into one big GMRF, on which inference can be carried out. Other features include functions to aid data import, md5 digest for
#' quick load following processing, standard pre-processing methods for observations, ability to cater for differing change of support, the possibility to
#' cater for "fine-scale variation", i.e. a fine-scale component over and above observation noise, domain boundaries for masking observations/inclusion
#' in observation processes, simultaneous combination of spatial/spatio-temporal fields, some functions for finding length scales from data for prior modelling and more....
#'
#' Note that if you want to use tif files for extracting some information e.g. topography, this library needs rgdal which may be awkward to install on Linux systems. First you need to install gdal, and then proj.4. I encountered some difficulties but the following links provided some answers.
#' http://trac.osgeo.org/proj/ticket/153
#' http://lightningismyname.blogspot.co.uk/2010/10/fixing-errors-with-nan-c-const-and-gcc.html
#' Note that if you are a local user you will need to alter the LD_LIBRARY_PATH, e.g.
#' export LD_LIBRARY_PATH="/home/glacio/ggazm/gdal/lib/:$HOME/proj/lib"
#' Also, you will need to point the install.packages function in the right direction e.g.
#' install.packages("rgdal",configure.args=c('--with-proj-include=$HOME/proj4/include','--with-proj-lib=$HOME/proj4/lib'))
#'
#'
#' @docType package
#' @useDynLib FRK
#' @import ggplot2
#' @import spam
#' @import Matrix
#' @import sp
#' @import parallel
#' @importClassesFrom gpclib gpc.poly
#' @importFrom raster raster
#' @importFrom gpclib area.poly
#' @importFrom akima interp
#' @importFrom digest digest
#' @importFrom maptools readShapeSpatial
#' @importFrom dplyr group_by summarise arrange
#' @importFrom plyr empty dlply . daply ddply adply rbind.fill
#' @importFrom deldir deldir
#' @importFrom SDMTools pnt.in.poly
#' @importFrom network network as.matrix.network.adjacency as.matrix.network.edgelist
#' @importFrom scales brewer_pal gradient_n_pal muted
#' @importFrom grid unit grid.newpage pushViewport viewport grid.layout
#' @importFrom geometry tsearch
#' @importFrom fields rdist fields.rdist.near
#' @importFrom actuar dinvgamma rinvgamma qinvgamma
#' @importFrom igraph graph.adjacency graph.bfs
#' @name FRK
NULL


