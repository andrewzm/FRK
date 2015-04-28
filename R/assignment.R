#' @title Assignment methods
#' @name [<--methods
#' @docType methods
#' @rdname assignment
#' @description Methods for \code{"[<-"}, i.e., extraction or subsetting of elements in the data frame of the MVST object
#' @examples 
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' Mesh["z"] <- sin(Mesh["x"]/1000)*cos(Mesh["y"]/1000)
NULL


setMethod("[",signature = "GMRF_basis",function(x,i,j) { return(x@G@rep[i][,])})

#' @rdname assignment
#' @aliases [<-,GMRF_basis,ANY,ANY-method
#' @export
setMethod("[<-",signature = "GMRF_basis",function(x,i,j,value) {
  x@G@rep[i] <- value
  return(x)
})
setMethod("[",signature = "GMRF",function(x,i,j) { return(x@rep[i][,])})

#' @rdname assignment
#' @aliases [<-,GMRF,ANY,ANY-method
#' @export
setMethod("[<-",signature = "GMRF",function(x,i,j,value) {
  x@rep[i] <- value
  return(x)
})
setMethod("[",signature = "Obs",function(x,i,j) { return(x@df[i][,])})

#' @rdname assignment
#' @aliases [<-,Obs,ANY,ANY-method
#' @export
setMethod("[<-",signature = "Obs",function(x,i,j,value) {
  x@df[i] <- value
  return(x)
})

setMethod(f="[", signature="Basis", definition=function(x,i,j) {return(x@pars$vars[i][,])})


#' @rdname assignment
#' @aliases [<-,Basis,ANY,ANY-method
#' @export
setMethod(f="[<-", signature="Basis", definition=function(x,i,j,value) {
  x@pars$vars[i] <- value
  return(x)})
