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

#' Global rhwrapper
#' @noRd
rhwrapper <- function(Ntot = 20, N = 10,type="data.frame",f_expr,...) {
    if(!(require("Rhipe"))) stop("Rhipe required for using Hadoop backend")

    njobs <- ceiling(Ntot/N)
    map1 <- expression({
        #suppressMessages(library(Matrix))
        #suppressMessages(library(sp))
        #suppressMessages(library(FRK))
        lapply(seq_along(map.keys), function(r) {
            #idx <- as.numeric((map.values[[r]]-1)*N + (1:N))
            #idx <- idx[which(idx <= Ntot)]
            outputkey <- map.values[[r]]
            #f <- eval(f_expr)
	    #outputvalue <- f(idx)
	    outputvalue <- data.frame(x = rnorm(10),y=rnorm(10))
            #outputvalue <- serialize(outputvalue,NULL)
            rhcollect(outputkey, outputvalue)
        })
    })

    mapreduce1 <- Rhipe::rhwatch(
        map      = map1,
        input    = njobs,
        output   = Rhipe::rhfmt("reduced_data", type = "sequence"),
        mapred = list(mapreduce.job.maps=njobs,
                      mapreduce.job.reduces=0),
        readback = TRUE#,
        #parameters = c(list(...),f_expr=f_expr,N=N,Ntot=Ntot,type=type)
    )

     Y <- do.call("rbind",lapply(mapreduce1,function(x) unserialize(x[[2]])))


    return(Y)
}

#' @title sp::over using Rhipe
#' @noRd
.rhover <- quote(function(idx) {
    suppressMessages(library(sp))
    sp::over(sp_pols[idx,],
              data_sp[c(av_var,"Nobs","std")],
              fn=sum)
})

#' @title eval_basis using Rhipe
#' @noRd
.rhpoint_eval_fn <- quote(function(idx) {
     suppressMessages(library(Matrix))
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

#' @title .dlply using Rhipe
#' @noRd
.rhdlply <- quote(function(idx) {
    suppressMessages(library(sp))
    plyr::dlply(df[idx,],keys,eval(dfun))
})
