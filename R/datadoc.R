#' @title AIRS data for May 2003
#'
#' @description Mid-troposheric CO2 measurements from the Atmospheric InfraRed Sounder (AIRS). The data represent measurements between 60 degrees S and 90 degrees N at roughly 1:30 pm local time on 1 May through 15 May 2003.
#'
#' @format A data frame with 209631 rows and 7 variables:
#' \describe{
#'   \item{year}{year of retrieval}
#'   \item{month}{month of retrieval}
#'   \item{day}{day of retrieval}
#'   \item{lat}{latitude coordinate of retrieval}
#'   \item{lon}{longitude coordinate of retrieval}
#'   \item{co2avgret}{CO2 mole fraction in ppm}
#'   \item{co2std}{standard deviation of CO2 reading in ppm}
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
#'   \item{lat}{latitude coordinate of retrieval}
#'   \item{lon}{longitude coordinate of retrieval}
#' }
#' @docType data
#' @references  Chahine, M. et al. (2006). AIRS: Improving weather forecasting and providing new data on greenhouse gases. Bulletin of the American Meteorological Society 87, 911--26.
"NOAA_df_1990"

#' @title ISEA Aperture 3 Hexagon (ISEA3H) Discrete Global Grid
#'
#' @description The data used here was obtained from http://webpages.sou.edu/~sahrk/dgg/isea.old/gen/isea3h.html and represent ISEA discrete global grids (DGGRIDs) generated using the \code{DGGRID} software. The original .gen files were converted to a data frame using the function \code{dggrid_gen_to_df}, available with the \code{FRK} package.
#' @format A data frame with 766280 rows and 5 variables:
#' \describe{
#'   \item{id}{grid identification number within the given resolution}
#'   \item{lon}{longitude coordinate}
#'   \item{lat}{latitude coordinate}
#'   \item{res}{DGGRID resolution (0 -- 9)}
#'   \item{centroid}{A 0-1 variable, indicating whether the point describes the centroid of the polygon, or is a point describing the polygon itself}
#' }
#' @docType data
#' @references  Sahr, K. (2008) Location coding on icosahedral aperture 3 hexagon discrete global grids. Computers, Environment and Urban Systems, 32(3):174-187.
"isea3h"
