#' @rdname nrow
#' @aliases nrow,GMRF-method
setMethod("nrow",signature("GMRF"),function(x){return (x@n)})

#' @rdname nrow
#' @aliases nrow,GMRF_basis-method
setMethod("nrow",signature("GMRF_basis"),function(x){return (nrow(x@G))})

#' @rdname nrow
#' @aliases nrow,FEBasis-method
setMethod("nrow",signature("FEBasis"),function(x){return (x@n)})

#' @rdname nrow
#' @aliases nrow,Obs-method
setMethod("nrow",signature("Obs"),function(x){return (nrow(x@df))})


#' @rdname subset
#' @aliases subset,Basis-method
setMethod("subset",signature="Basis",definition=function(x,...) {
  return(subset(x@pars$vars,...))
})

#' @rdname subset
#' @aliases subset,GMRF_basis-method
setMethod("subset",signature="GMRF_basis",definition=function(x,...) {
  return(subset(x@G,...))
})

#' @rdname subset
#' @aliases subset,GMRF-method
setMethod("subset",signature="GMRF",definition=function(x,...) {
  return(subset(x@rep,...))
})

#' @rdname subset
#' @aliases subset,Obs-method
setMethod("subset",signature = "Obs",function(x,...) {
  x@df <- subset(x@df,...)
  return(x)
})

#' @rdname subset
#' @aliases subset,Obs_poly-method
setMethod("subset",signature = "Obs_poly",function(x,...) {
  x@df <- subset(x@df,...)
  x@df2 <- subset(x@df2,...)
  x@pol <- x@pol[x@df$n]
  if ("P" %in% names(x@args)) x <- setalpha(x,x@args$alpha0,x@args$av_dist)
  return(x)
})

#' @rdname head
#' @aliases head,block-method
setMethod("head",signature="block",function(x,...) { return(head(getDf(x),...))})

#' @rdname head
#' @aliases head,Basis-method
setMethod("head",signature="Basis",function(x,...) { return(head(getDf(x),...))})

#' @rdname tail
#' @aliases tail,block-method
setMethod("tail",signature="block",function(x,...) { return(tail(getDf(x),...))})

#' @rdname tail
#' @aliases tail,Basis-method
setMethod("tail",signature="Basis",function(x,...) { return(tail(getDf(x),...))})

