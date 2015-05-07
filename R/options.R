defaults = list(progress = TRUE, verbose = FALSE, parallel=6L)

opts_FRK <- list(
    set = function(opt,value) {
        defaults[opt] <<- value
    },
    get = function(opt) {
        defaults[opt]
    }
)
