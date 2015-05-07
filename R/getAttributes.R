#' @rdname nbasis
#' @aliases nbasis,Basis-method
setMethod("nbasis",signature(.Object="Basis"),function(.Object) {return(.Object@n)})

#' @rdname nbasis
#' @aliases nbasis,SRE-method
setMethod("nbasis",signature(.Object="SRE"),function(.Object) {return(nbasis(.Object@basis))})

#' @rdname type
#' @aliases type,manifold-method
setMethod("type",signature(.Object="manifold"),function(.Object) {
    return(.Object@type)
})


#' @rdname manifold
#' @aliases manifold,Basis-method
setMethod("manifold",signature(.Object="Basis"),function(.Object) {
    return(.Object@manifold)
})


setMethod("getData",signature(.Object="Obs"),function(.Object) {
  return(.Object@df$z)
})

#' @rdname getDf
#' @aliases getDf,Graph-method
setMethod("getDf",signature(.Object="Graph"),function(.Object) {
  Olist <- extractClass(.Object@v,"Obs")
  if(length(Olist)>1) stop("Only works with one observation for now")
  return(Olist[[1]]@df)
})

#' @rdname getDf
#' @aliases getDf,FEBasis-method
setMethod("getDf",signature(.Object="FEBasis"),function(.Object) {
  return(.Object@pars$vars)
})

#' @rdname getDf
#' @aliases getDf,GMRF-method
setMethod("getDf",signature(.Object="GMRF"),function(.Object) {
  return(.Object@rep)
})

#' @rdname getDf
#' @aliases getDf,GMRF_basis-method
setMethod("getDf",signature(.Object="GMRF_basis"),function(.Object) {
  return(.Object@G@rep)
})

#' @rdname getDf
#' @aliases getDf,Obs-method
setMethod("getDf",signature(.Object="Obs"),function(.Object) {
  return(.Object@df)
})

#' @rdname getPrecision
#' @aliases getPrecision,Obs-method
setMethod("getPrecision",signature(.Object="Obs"),function(.Object) {
  return(sparsediag(1/.Object@df$std^2))
})

#' @rdname getPrecision
#' @aliases getPrecision,GMRF-method
setMethod("getPrecision",signature(.Object="GMRF"),function(.Object) {
  return(.Object@Q)
})

#' @rdname getPrecision
#' @aliases getPrecision,GMRF_basis-method
setMethod("getPrecision",signature(.Object="GMRF_basis"),function(.Object) {
  return(getPrecision(.Object@G))
})

#' @rdname getMean
#' @aliases getMean,GMRF-method
setMethod("getMean",signature(.Object="GMRF"),function(.Object) {
  return(.Object@mu)
})

#' @rdname getMean
#' @aliases getMean,GMRF_basis-method
setMethod("getMean",signature(.Object="GMRF_basis"),function(.Object) {
  return(getMean(.Object@G))
})

#' @rdname getC
#' @aliases getC,Graph_2nodes-method
setMethod("getC",signature(.Object="Graph_2nodes"),function(.Object) {
 return(getC(.Object@e[[1]]))
})


#' @rdname getC
#' @aliases getC,linkGO-method
setMethod("getC",signature(.Object="linkGO"),function(.Object) {
  return(.Object@Cmat)
})
