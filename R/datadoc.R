#' ICESAT data.
#' 
#' A sample of the Ice, Cloud and Elevation dataset, where each point represents an average trend in m/yr over a 20km x 20km grid box. The data has the following fields
#' 
#' \itemize{
#'   \item t. year of observation (0--6 where 0 is 2003, and 6 is 2009)
#'   \item x. x coordinate in km (under polar stereographic projection)
#'   \item y. y coordinate in km (under polar stereographic projection)
#'   \item std. std standard error of trend
#'   \item z. trend
#'   \item obs_name. the name of the observation
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 10000 rows and 6 variables
#' @name icesat
NULL

#' Example mesh
#' 
#' A list containing the mesh details of the Antarctic ice sheet at moderate (~200km) resolution
#' 
#' \itemize{
#'   \item p. a 1984 x 2 matrix containing the vertices
#'   \item t. a 3911 x 3 matrix containing the trianulations
#'   \item M. a 1984 x 1984 sparse mass matrix
#'   \item K. a 1984 x 1984 sparse stiffness matrix
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A list with four fields
#' @name surf_fe
NULL

#' Boundaries for Antarctica
#' 
#' A list of data frames containing useful shapefiles for plotting results relating to Antarctica.
#' 
#' \itemize{
#'   \item grounding_sub. a 1984 x 2 matrix containing the vertices
#'   \item coast_sub. a 3911 x 3 matrix containing the trianulations
#'   \item Islands. a 1984 x 1984 sparse mass matrix
#'   \item K. a 1984 x 1984 sparse stiffness matrix
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A list with three fields. Each data frame has at least the following fields
#'  \itemize{
#'   \item x. x coordinate in km under a polar stereographic projection
#'   \item y. y coordinate in km under a polar stereographic projection
#'   \item id. polygon id
#' }
#' @name shapefiles
NULL