#' @title AIRS data for May 2003
#'
#' @description Mid-tropospheric CO2 measurements from the Atmospheric InfraRed Sounder (AIRS). The data are measurements between 60 degrees S and 90 degrees N at roughly 1:30 pm local time on 1 May through to 15 May 2003.
#'
#' @format A data frame with 209631 rows and 7 variables:
#' \describe{
#'   \item{year}{year of retrieval}
#'   \item{month}{month of retrieval}
#'   \item{day}{day of retrieval}
#'   \item{lon}{longitude coordinate of retrieval}
#'   \item{lat}{latitude coordinate of retrieval}
#'   \item{co2avgret}{CO2 mole fraction retrieval in ppm}
#'   \item{co2std}{measurement error of CO2 retrieval in ppm}
#' }
#' @docType data
#' @references  Chahine, M. et al. (2006). AIRS: Improving weather forecasting and providing new data on greenhouse gases. Bulletin of the American Meteorological Society 87, 911--26.
"AIRS_05_2003"

#' @title NOAA maximum temperature data for 1990-1993
#'
#' @description Maximum temperature data obtained from the National Oceanic and Atmospheric Administration (NOAA) for a part of the USA between 1990 and 1993 (inclusive).
#'
#' @format A data frame with 196253 rows and 8 variables:
#' \describe{
#'   \item{year}{year of retrieval}
#'   \item{month}{month of retrieval}
#'   \item{day}{day of retrieval}
#'   \item{z}{dependent variable}
#'   \item{proc}{variable name (Tmax)}
#'   \item{id}{station id}
#'   \item{lon}{longitude coordinate of retrieval}
#'   \item{lat}{latitude coordinate of retrieval}
#' }
#' @docType data
#' @references  Chahine, M. et al. (2006). AIRS: Improving weather forecasting and providing new data on greenhouse gases. Bulletin of the American Meteorological Society 87, 911--26.
"NOAA_df_1990"

#' @title ISEA Aperture 3 Hexagon (ISEA3H) Discrete Global Grid
#'
#' @description The data used here was obtained from http://webpages.sou.edu/~sahrk/dgg/isea.old/gen/isea3h.html and represent ISEA discrete global grids (DGGRIDs) generated using the \code{DGGRID} software. The original .gen files were converted to a data frame using the function \code{dggrid_gen_to_df}, available with the \code{dggrids} package. Only resolutions 0--6 are supplied with \code{FRK} and note that resolution 0 of ISEA3H is equal to resolution 1 in \code{FRK}. For higher resolutions please install \code{dggrids} from \code{https://github.com/andrewzm/dggrids}.
#' @format A data frame with 284208 rows and 5 variables:
#' \describe{
#'   \item{id}{grid identification number within the given resolution}
#'   \item{lon}{longitude coordinate}
#'   \item{lat}{latitude coordinate}
#'   \item{res}{DGGRID resolution (0 -- 6)}
#'   \item{centroid}{A 0-1 variable, indicating whether the point describes the centroid of the polygon, or whether it is a boundary point of the polygon}
#' }
#' @docType data
#' @references  Sahr, K. (2008) Location coding on icosahedral aperture 3 hexagon discrete global grids. Computers, Environment and Urban Systems, 32, 174--187.
"isea3h"

#' @title World map
#'
#' @description This world map was extracted from the package \code{maps} v.3.0.1 by running \code{map_data("world")}. To reduce the data size, only every third point of this data frame is contained in \code{worldmap}.
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
#' @references  Original S code by Becker, R.A. and Wilks, R.A. This R version is by Brownrigg, R. Enhancements have been made by Minka, T.P. and Deckmyn, A. (2015) maps: Draw Geographical Maps, R package version 3.0.1.
"worldmap"
