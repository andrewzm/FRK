#' @title Options for FRK package
#' @export
opts_FRK <- list(
    set = function(opt,value) {
        if(!(opt %in% c("progress","verbose","parallel"))) stop("opt needs to be one of ('progress','verbose','parallel')")
        defaults[opt] <<- value
    },
    get = function(opt) {
        if(!(opt %in% c("progress","verbose","parallel"))) stop("opt needs to be one of ('progress','verbose','parallel')")
        defaults[opt]
    }
)


defaults = list(progress = TRUE, verbose = FALSE, parallel=6L)


