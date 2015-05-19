#' @title Options for FRK package
#' @export
new_opts_FRK <- function(d = list(progress = TRUE, verbose = FALSE, parallel=6L, Rhipe=FALSE)) {
    defaults = d
    list(set = function(opt,value) {
        if(!(opt %in% c("progress","verbose","parallel","Rhipe"))) stop("opt needs to be one of ('progress','verbose','parallel')")
        defaults[[opt]] <<- value
    },
    get = function(opt) {
        if(!(opt %in% c("progress","verbose","parallel","Rhipe"))) stop("opt needs to be one of ('progress','verbose','parallel')")
        defaults[[opt]]
    }
)}

opts_FRK = new_opts_FRK()






