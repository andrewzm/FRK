## initialise options: Run when FRK is loaded
new_opts_FRK <- function(d = list(progress = TRUE, verbose = FALSE, parallel=1L)) {
    defaults = d                               # detauls to showing progress bars, no verbose and single core computation
    ## there are two functions in this list: $set() and $get()
    list(set = function(opt = c("progress","verbose","parallel"),value) {
        opt <- match.arg(opt)                 # check option being set
        value <- .option_check(opt,value)     # check value is appropriate
        defaults[[opt]] <<- value             # set value

        if(opt == "parallel") {               # if switching off or on a cluster, take the appropriate steps
            if(!is.null(defaults[["cl"]])) {
                stopCluster(defaults[["cl"]])
                defaults[["cl"]] <<- NULL     # we are updating defaults in the parent environment
            }
            if (value > 1) {
                defaults[["cl"]] <<- makeCluster(value,useXDR=FALSE)
            }
        }
    },
    ## Second function (get)
    get = function(opt = c("progress","verbose","parallel","cl")) {
        opt <- match.arg(opt)
        defaults[[opt]]  # just retrieve option
    })
}


#' @title FRK options
#' @description The main options list for the FRK package.
#' @format List of 2
#' \itemize{
#'   \item{\code{$}  }{\code{set:function(opt,value)}}
#'   \item{\code{$}  }{\code{get:function(opt)}}
#' }
#' @details \code{opts_FRK} is a list containing two functions, \code{set} and \code{get}, which can be used to set options and retrieve options, respectively. Currently \code{FRK} uses three options:
#' \itemize{
#'  \item{"progress":}{ a flag indicating whether progress bars should be displayed or not}
#'  \item{"verbose":}{ a flag indicating whether certain progress messages should be shown or not}
#'  \item{"parallel":}{ an integer indicating the number of cores to use. A number 0 or 1 indicates no parallelism}
#' }
#' @examples
#' opts_FRK$set("progress",1L)
#' opts_FRK$get("parallel")
#' @export
opts_FRK = new_opts_FRK()

#####################
## NOT EXPORTED #####
#####################

## Checking of options. The code is self-explanatory
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
    }

    if(opt == "parallel")
        if(!(value >= 0))
            stop("parallel should be a nonnegative integer")
    if(opt == "parallel")
        if(!requireNamespace("parallel"))
            stop("package parallel is required for using multiple cores. Please install parallel")

    value

}






