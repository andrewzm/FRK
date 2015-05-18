#' Global rhwrapper
#' @noRd
rhwrapper <- function(Ntot = 20, N = 10,type="data.frame",f_expr,...) {
    njobs <- ceiling(Ntot/N)

    map1 <- expression({
        suppressMessages(library(Matrix))
        suppressMessages(library(sp))
        library(FRK)
        lapply(seq_along(map.keys), function(r) {
            idx <- as.numeric((map.values[[r]]-1)*N + (1:N))
            idx <- idx[which(idx <= Ntot)]
            outputkey <- map.values[[r]]
            f <- eval(f_expr)
            outputvalue <- f(idx)
            outputvalue <- serialize(outputvalue,NULL)
            rhcollect(outputkey, outputvalue)
        })
    })

    mapreduce1 <- rhwatch(
        map      = map1,
        input    = njobs,
        output   = rhfmt("reduced_data", type = "sequence"),
        mapred   = list(mapred.map.tasks=njobs,
			mapred.reduce.tasks=0),
        readback = TRUE,
        parameters = c(list(...),f_expr=f_expr,N=N,Ntot=Ntot,type=type)
    )

    if(type=="Matrix") {
        Y <- do.call("rBind",lapply(mapreduce1,function(x) unserialize(x[[2]])))
     } else {
        Y <- do.call("rbind",lapply(mapreduce1,function(x) unserialize(x[[2]])))
     }

    Y
}

#' @title sp::over using Rhipe
#' @noRd
.rhover <- quote(function(idx) {
    sp::over(sp_pols[idx,],
              data_sp[c(av_var,"Nobs","std")],
              fn=sum)
})

#' @title eval_basis using Rhipe
#' @noRd
.rhpoint_eval_fn <- quote(function(idx) {
    ## x <- sapply(flist,function(fn) fn(s[idx,,drop=FALSE]))
     x <- sapply(seq_along(flist), function(i) {
                e <- envlist[[i]]
               c <- e$c
                manifold <- e$manifold
                R <- e$R
                fn <- eval(flist[[i]])
                stopifnot(is(fn,"function"))
                fn(s[idx,,drop=FALSE])
     })
    Matrix(as(x,"matrix"))
})

#' @title SRE.predict using Rhipe
#' @noRd
.rhSRE.predict <- quote(function(idx) {
   FRK:::.SRE.predict(Sm=Sm,pred_locs=pred_locs[idx,],use_centroid=use_centroid)
})

