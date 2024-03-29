% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGeneric.R, R/basisfns.R
\name{nres}
\alias{nres}
\alias{nres,Basis-method}
\alias{nres_basis,Basis-method}
\alias{nres,TensorP_Basis-method}
\alias{nres,SRE-method}
\alias{nres_SRE,SRE-method}
\title{Return  the number of resolutions}
\usage{
nres(b)

\S4method{nres}{Basis}(b)

\S4method{nres}{TensorP_Basis}(b)

\S4method{nres}{SRE}(b)
}
\arguments{
\item{b}{object of class \code{Basis} or \code{SRE}}
}
\description{
Return the number of resolutions from a basis function object.
}
\examples{
library(sp)
set.seed(1)
d <- data.frame(lon = runif(n=500,min = -179, max = 179),
                lat = runif(n=500,min = -90, max = 90),
                z = rnorm(500))
coordinates(d) <- ~lon + lat
slot(d, "proj4string") = CRS("+proj=longlat")

### Now create basis functions on sphere
G <- auto_basis(manifold = sphere(),data=d,
                nres = 2,prune=15,
                type = "bisquare",
                subsamp = 20000)
nres(G)
}
\seealso{
\code{\link{auto_basis}} for automatically constructing basis functions and \code{\link{show_basis}} for visualising basis functions.
}
