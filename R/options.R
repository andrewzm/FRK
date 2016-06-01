## initialise options: Run when FRK is loaded
new_opts_FRK <- function(d = list(progress = TRUE, verbose = FALSE, parallel=1L)) {
    defaults = d
    list(set = function(opt,value) {
        if(!(opt %in% c("progress","verbose","parallel")))
            stop("opt needs to be one of ('progress','verbose','parallel')")
        value <- .option_check(opt,value)
        defaults[[opt]] <<- value
    },
    get = function(opt) {
        if(!(opt %in% c("progress","verbose","parallel","Rhipe")))
            stop("opt needs to be one of ('progress','verbose','parallel')")
        defaults[[opt]]
    }
    )}

.option_check <- function(opt,value) {
    if(opt == "progress")
        if(!(value == TRUE | value == FALSE))
            stop("progress should be TRUE or FALSE")

    if(opt == "verbose")
        if(!(value == TRUE | value == FALSE))
            stop("verbose should be TRUE or FALSE")

    if(opt == "parallel") {
        if(!is.integer(value))
            stop("parallel should be a nonnegative integer")
        if(Sys.info()[['sysname']] == "Windows") {
            if(!value == 1L) {
                warning("Windows detected. Currently FRK is only parallelisable on
                        Linux/Mac systems. Coercing number of cores to 1")
                value <- 1L
            }
        }
    }

    if(opt == "parallel")
        if(!(value >= 0))
            stop("parallel should be a nonnegative integer")
    if(opt == "parallel")

    if(!requireNamespace("parallel"))
        stop("package parallel is required for using multiple cores. Please install parallel")

    value

}

#' @title FRK options
#' @description The main options list for the FRK package.
#' @format List of 2
#' \itemize{
#'   \item{\code{$}  }{\code{set:function(opt,value)}}
#'   \item{\code{$}  }{\code{get:function(opt)}}
#' }
#' @details \code{opts_FRK} is a list containing two functions, \code{set} and \code{get}, which can be used to set options and retrieve options, respectively. Currently \code{FRK} uses four options:
#' \itemize{
#'  \item{"progress":}{ a flag indicating whether progress bars should be displayed or not}
#'  \item{"verbose":}{ a flag indicating whether certain progress monitors should be shown or not}
#'  \item{"parallel":}{ an integer indicating the number of cores to use. A number 0 or 1 indicates no parallelism}
#' }
#' @examples
#' opts_FRK$set("parallel",2L)
#' opts_FRK$get("parallel")
#' @export
opts_FRK = new_opts_FRK()






