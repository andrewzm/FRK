# FRK: An R Software package for spatial and spatio-temporal prediction
# with large datasets.
# Copyright (c) 2017 University of Wollongong
# Author: Andrew Zammit-Mangion, azm (at) uow.edu.au
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#' @title AIRS data for May 2003
#' @description Mid-tropospheric CO2 measurements from the Atmospheric InfraRed Sounder (AIRS).
#' The data are measurements between 60 degrees S and 90 degrees N at roughly 1:30 pm local
#' time on 1 May through to 15 May 2003. (AIRS does not release data below 60 degrees S.)
#' @format A data frame with 209631 rows and 7 variables:
#' \describe{
#'   \item{year}{year of retrieval}
#'   \item{month}{month of retrieval}
#'   \item{day}{day of retrieval}
#'   \item{lon}{longitude coordinate of retrieval}
#'   \item{lat}{latitude coordinate of retrieval}
#'   \item{co2avgret}{CO2 mole fraction retrieval in ppm}
#'   \item{co2std}{standard error of CO2 retrieval in ppm}
#' }
#' @docType data
#' @references  Chahine, M. et al. (2006). AIRS: Improving weather forecasting and
#' providing new data on greenhouse gases. Bulletin of the American Meteorological
#' Society 87, 911--26.
"AIRS_05_2003"

#' @title NOAA maximum temperature data for 1990--1993
#' @description Maximum temperature data obtained from the National Oceanic and Atmospheric
#' Administration (NOAA) for a part of the USA between 1990 and 1993 (inclusive).
#' See https://iridl.ldeo.columbia.edu/ SOURCES/.NOAA/.NCDC/.DAILY/.FSOD/.
#' @format A data frame with 196,253 rows and 8 variables:
#' \describe{
#'   \item{year}{year of retrieval}
#'   \item{month}{month of retrieval}
#'   \item{day}{day of retrieval}
#'   \item{z}{dependent variable}
#'   \item{proc}{variable name (Tmax)}
#'   \item{id}{station id}
#'   \item{lon}{longitude coordinate of measurement station}
#'   \item{lat}{latitude coordinate of measurement station}
#' }
#' @docType data
#' @references  National Climatic Data Center, March 1993: Local Climatological Data.
#' Environmental Information summary (C-2), NOAA-NCDC, Asheville, NC.
"NOAA_df_1990"

#' @title ISEA Aperture 3 Hexagon (ISEA3H) Discrete Global Grid
#'
#' @description The data used here were obtained from
#' https://webpages.sou.edu/~sahrk/dgg/isea.old/gen/isea3h.html and represent ISEA
#' discrete global grids (DGGRIDs) generated using the \code{DGGRID} software.
#' The original .gen files were converted to a data frame using the function \code{dggrid_gen_to_df},
#' available with the \code{dggrids} package. Only resolutions 0--6 are supplied with \code{FRK}
#' and note that resolution 0 of ISEA3H is equal to resolution 1 in \code{FRK}. For higher
#' resolutions \code{dggrids} can be installed from \code{https://github.com/andrewzm/dggrids/}
#' using \code{devtools}.
#' @format A data frame with 284,208 rows and 5 variables:
#' \describe{
#'   \item{id}{grid identification number within the given resolution}
#'   \item{lon}{longitude coordinate}
#'   \item{lat}{latitude coordinate}
#'   \item{res}{DGGRID resolution (0 -- 6)}
#'   \item{centroid}{A 0-1 variable, indicating whether the point describes the centroid of the polygon,
#'   or whether it is a boundary point of the polygon}
#' }
#' @docType data
#' @references  Sahr, K. (2008). Location coding on icosahedral aperture 3 hexagon discrete global grids. Computers, Environment and Urban Systems, 32, 174--187.
"isea3h"

#' @title World map
#' @description This world map was extracted from the package \code{maps} v.3.0.1 by
#' running \code{ggplot2::map_data("world")}. To reduce the data size, only every third point of
#' this data frame is contained in \code{worldmap}.
#' @format A data frame with 33971 rows and 6 variables:
#' \describe{
#'   \item{long}{longitude coordinate}
#'   \item{lat}{latitude coordinate}
#'   \item{group}{polygon (region) number}
#'   \item{order}{order of point in polygon boundary}
#'   \item{region}{region name}
#'   \item{subregion}{subregion name}
#' }
#' @docType data
#' @references  Original S code by Becker, R.A. and Wilks, R.A. This R version is by
#' Brownrigg, R. Enhancements have been made by Minka, T.P. and Deckmyn, A. (2015)
#' maps: Draw Geographical Maps, R package version 3.0.1.
"worldmap"



#' @title MODIS cloud data
#' @description An image of a cloud taken by the Moderate Resolution
#' Imaging Spectroradiometer (MODIS) instrument aboard the Aqua satellite (MODIS
#' Characterization Support Team, 2015). 
#' @format A data frame with 33,750 rows and 3 variables:
#' \describe{
#'   \item{x}{x-coordinate}
#'   \item{y}{y-coordinate}
#'   \item{z}{binary dependent variable: 1 if cloud is present, 0 if no cloud. This
#'   variable has been thresholded from the original continuous measurement of
#'   radiance supplied by the MODIS instrument}
#'   \item{z_unthresholded}{The original continuous measurement of
#'   radiance supplied by the MODIS instrument}
#' }
#' @docType data
#' @references  MODIS Characterization Support Team (2015). MODIS 500m Calibrated Radiance Product.NASA MODIS Adaptive Processing System, Goddard Space Flight Center, USA.
"MODIS_cloud_df"


#' @title Americium soil data
#' @description Americium (Am) concentrations in a spatial domain immediately surrounding the location at which nuclear devices were detonated at Area 13 of the Nevada Test Site, between 1954 and 1963. 
#' @format A data frame with 212 rows and 3 variables:
#' \describe{
#'   \item{Easting}{Easting in metres}
#'   \item{Northing}{Northing in metres}
#'   \item{Am}{Americium concentration in 1000 counts per minute}
#' }
#' @docType data
#' @references Paul R, Cressie N (2011). “Lognormal block kriging for contaminated soil.” European Journal of Soil Science, 62, 337–345.
"Am_data"
