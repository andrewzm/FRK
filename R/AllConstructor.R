#### CONSTRUCTORS ###########

#' @title manifold
#'
#' @description This function initialises a manifold
setMethod("initialize",signature="manifold",function(.Object) {
    ## General manifold checks can come in here
    .Object
})


#' @title sphere
#'
#' @description This function initialises a sphere

sphere <- function(radius=1,metric=gc_dist(R=radius)) {
  stopifnot(dimensions(metric)==2L)
  stopifnot(radius>0)
  new("sphere",metric=metric,radius=radius)
}

setMethod("initialize",signature="sphere",function(.Object,radius=1,metric=gc_dist(R=radius)) {
    .Object@type <- "sphere"
    .Object@metric <- metric
    .Object@radius <- radius
    callNextMethod(.Object)})

#' @title plane
#'
#' @description This function initialises a plane
plane <- function(metric=Euclid_dist(dim=2L)) {
    stopifnot(dimensions(metric)==2L)
    new("plane",metric=metric)
}

setMethod("initialize",signature="plane",function(.Object,metric=Euclid_dist(dim=2L)) {
    .Object@type <- "plane"
    .Object@metric <- metric
    callNextMethod(.Object)})



#' @title real_line
#'
#' @description This function initialises a sphere
real_line <- function(metric=Euclid_dist(dim=1L)) {
    stopifnot(dimensions(metric)==1L)
    new("real_line",metric=metric)
}

setMethod("initialize",signature="real_line",function(.Object,metric=Euclid_dist(dim=1L)) {
    .Object@type <- "real_line"
    .Object@metric <- metric
    callNextMethod(.Object)})



#' @title Great-circle distance
#'
#' @description This function initialises a great-circle distance measure
gc_dist <- function(R=NULL) {
    new("measure",dist=function(x1,x2)  fields::rdist.earth(x1,x2,miles=F,R=R),dim=2L)
}

#' @title Euclidean distance
#'
#' @description This function initialises a Euclidean distance measure on 2D
Euclid_dist <- function(dim=2L) {
    stopifnot(is.integer(dim))
    new("measure",dist=function(x1,x2)  fields::rdist(x1,x2), dim=dim)
}

#' @title domain
#'
#' @description This function initialises a boundary polygon
domain <- function(m=sphere(),bndary=matrix(c(0,1,0,1),2,2)) {
    stopifnot(nrow(bndary) >= 2)
    stopifnot(ncol(bndary) == dimensions(m))
    new("domain",m,bndary)
}

setMethod("initialize",signature="domain",function(.Object,m = sphere(), bndary=matrix(c(0,1,0,1),2,2)) {
    .Object@m <- m
    .Object@bndary <- bndary
    .Object
    })

SRE <- function(f,data,basis) {
   L <- gstat:::gstat.formula(f,data=data)
   new("SRE",
       data=data,
       basis=basis,
       S = eval_basis(basis, s = coordinates(data)),
       V = Diagonal(x=data$std^2),
       Z = Matrix(L$y),
       X = as(L$X,"Matrix"))
}


#' @title GMRF
#'
#' @description This function initialises a GMRF with some mean \code{mu} and precision matrix \code{Q}. The returned object is of class \code{GMRF}
#'
#' @param mu mean, of class \code{matrix}
#' @param Q sparse precision matrix (of class \code{Matrix})
#' @param intrinsic set intrinsic level if Q is singular
#' @param n number of nodes. Note that this can be different from \code{nrow(mu)} if the system is multi-variate
#' @param name name of GMRF
#' @param rep data frame of length \code{N} with more details (for example axis, covariate information)
#' @param t_axis this is the time horizon for the system under consideration. If you are considering a spatial problem set this to zero.
#' @return Object of class GMRF
#' @keywords GMRF
#' @export
#' @examples
#'
#' require(Matrix)
#' # Create a GMRF
#' Q <- sparseMatrix(i=c(1,2,1,2),j=c(1,1,2,2),x=c(1,0.1,0.1,1))
#' mu <- matrix(c(0,0))
#' my_GMRF <- GMRF(mu=mu, Q=Q,name="my_first_GMRF")
#' print(getPrecision(my_GMRF))
#' print(getMean(my_GMRF))
#' print(getDf(my_GMRF))
GMRF <- function(mu=NA,Q= sparseMatrix(i=c(1,2),j=c(1,2),x=4),intrinsic=0,n=NULL,t_axis=0,
                 rep=data.frame(),name="none") {

  if(any(is.na(mu))) {
    mu <- matrix(0,nrow(Q),1)
  }
  stopifnot(nrow(mu) == nrow(Q))
  return(new("GMRF",mu=mu,Q=Q,intrinsic=intrinsic,n=n,t_axis=t_axis,rep=rep,name=name))
}

#' @title Random walk GMRF
#'
#' @description This function initialises a random walk and represents it as a Gaussian Markov Random Field with mean \code{mu} and precision matrix \code{Q}. Only random walks along the real line, first-order and second-order variants are implemented for now. As with other GMRFs, the user can specify a name. Also a data frame can be specified with more details on the GMRF.
#'
#' @param n number of vertices
#' @param order 1 or 2, depending on the order of the random walk
#' @param precinc precision constant (multiples the template precision matrix)
#' @param df data frame of length \code{n} with more details (for example axis, covariate information)
#' @param name name of GMRF
#' @return Object of class GMRF with zero mean
#' @keywords GMRF, random walk
#' @export
#' @examples
#'
#' require(Matrix)
#' # Createa a first-order random walk GMRF
#' my_RW <- GMRF_RW(n=10, order=1, precinc =2, name="my_first_RW")
#' print(getPrecision(my_RW))
#' print(getMean(my_RW))
#' print(getDf(my_RW))
GMRF_RW <- function(n = 10,order=1,precinc = 1,df=data.frame(),name="none") {

  stopifnot(order %in% c(1,2))

  if(is.null(n)) n<- nrow(mu)

  mu <- matrix(0,nrow=n,ncol=1)
  i = c(1:(n-1),1:(n-1))
  j = c(1:(n-1),2:n)
  x <- numeric(length=((n-1)*2))
  x[1:(n-1)] = -1
  x[n:((2*n)-2)] = 1
  Dmat = sparseMatrix(i,j,x=x)
  R = t(Dmat)%*%Dmat
  if (order == 1) {
    Q = precinc*R
    intrinsic = 1
  }
  if (order == 2) {
    R <- R %*% R
    R[1,(1:3)] = c(1,-2,1)
    R[2,(1:4)] = c(-2,5,-4,1)
    R[(n-1),(n-3):n] = c(1,-4,5,-2)
    R[(n),(n-2):n] = c(1,-2,1)
    Q <- precinc*R
    intrinsic = 2
  }

  return(new("GMRF",mu,Q,intrinsic=1,n=n,t_axis=0:(n-1),rep=df,name=name))
}

#' @title Vector auto-regressive model with GMRF represenation
#'
#' @description This function initialises a vector auto-regressive model and represents it as a Gaussian Markov Random Field with mean \code{mu} and precision matrix \code{Q}. This constructor differs from other GMRF constructors in that it takes function inputs
#' to define temporally evolving characteristics. The default representation is \deqn{x_{k+1} = \mu_k + A_kx_k + B_k\beta_k + e_k} where
#' \eqn{e_k \sim \mathcal{N}(0,Q_k)}. Note that in addition to covariates, a known mean \eqn{\mu_k} can be added, this can be omitted and
#' replaced appropriately with entries in \eqn{B_k}. A multi-variate vector auto-regressive model can be speficied by letting
#' \code{A_fun} and \code{Qw_fun} return matrices a multiple of the dimension of the underlying basis over which the GMRF is defined.
#'
#' @param mu_fun function of time \eqn{k}, returns matrix of size \eqn{n}
#' @param A_fun function of time \eqn{k}, returns sparse matrix of size \eqn{n\times n}
#' @param B_fun function of time \eqn{k}, returns  sparse matrix of size \eqn{n\times m}
#' @param Qw_fun function of time \eqn{k}, returns  sparse matrix of size \eqn{n\times n}
#' @param Qb prior precision matrix of \eqn{\beta}; sparse matrix of size \eqn{m \times m}
#' @param t_axis time axis of process
#' @param name name of VAR
#' @return Object of class VAR_Gauss which inherits from class \code{GMRF}.
#' @keywords auto-regressive model, multi-variate model
#' @export
#' @examples
#'
#' require(Matrix)
#' t_axis <- 0:10
#' mu <- function(k) return(matrix(0,length(t_axis),1))
#' A <- function(k)  return(sparsediag(0.4))
#' B <- function(k)  cBind(Imat(1),k*Imat(1))
#' Q <- function(k)  return(sparsediag(1))
#' Qb = bdiag(Imat(1),Imat(1))
#' VAR <- VAR_Gauss( mu_fun = mu,A=A, B=B, Qw = Q,t_axis = t_axis,Qb=Qb,name="firstVAR")
VAR_Gauss <- function(mu_fun = function(k) return(matrix(0,2,5)),
                      A_fun = function(k) return(Imat(2)),
                      B_fun =  function(k) return(emptySp()),
                      Qw_fun = function(k) return(Imat(2)),
                      t_axis = c(0:6),
                      Qb = emptySp(),
                      name="none")  {

  return( new("VAR_Gauss",mu_fun=mu_fun,A_fun=A_fun,B_fun=B_fun,Qw_fun=Qw_fun,t_axis=t_axis,Qb=Qb,name=name))

}



#' @title GMRF function over basis
#'
#' @description This function initialises an object of class \code{GMRF_basis} which defines a GMRF over a set of basis functions.
#'
#' @param G an object of class \code{GMRF}
#' @param Basis an object of class \code{Basis}
#' @return Object of class \code{GMRF_basis} (which inherits from class \code{process} and is thus also a process block)
#' @keywords GMRF, basis functions
#' @export
#' @examples
#'
#' G <- GMRF_RW(n=9)
#' Basis <-  initGRBFbasis(x = c(0,1), y = c(0,1), std=0.1,nx=9,ny=1)
#' print(GMRF_basis(G,Basis))
GMRF_basis <- function(G=new("GMRF"),Basis=new("Basis",pars=list(vars=data.frame(x=c(1,2))))) {
  stopifnot(is(G,"GMRF"))
  stopifnot(is(Basis,"Basis"))
  return(new("GMRF_basis",G=G,Basis=Basis))
}

#' @title Observation block
#'
#' @description This function initialises an object of class \code{Obs} which defines a an observation data set. By default, this is for observations with negligible spatial footprint. For larger supports, use \code{Obs_poly}.
#'
#' @param df a data frame which should contain at least 5 entries, \code{x,y,t,z} and \code{std} which denote the horizontal, vertical and temporal indices of the observations, the value and error respectively. Alternatively this could be a path name.
#' @param name the name of the observation process
#' @param remove_cross_ins removes data outside a circle centred at zero with specified radius. Convenient when working with satellite data in polar stereographic projection when some cross-ins are detected.
#' @param ... other arguments passed on to \code{preprocess_obs}
#' @return Object of class \code{Obs} (which inherits from class \code{block} and is thus also a block)
#' @keywords Observations, change of support, block
#' @export
#' @examples
#' O <- Obs(df=data.frame(x=runif(5),y=runif(5),t=c(1,1,1,2,2),z=runif(5),std=runif(5)))
#' print(O)
#' plot(subset(O,t==1),"z",pt_size=4)
Obs <- function(df,name="Obs",remove_cross_ins=0,...) {
    return(new("Obs",df=df,name=name,remove_cross_ins=remove_cross_ins,...))
}

#' @title Observation block with support
#'
#' @description This function initialises an object of class \code{Obs_poly} which defines a an observation data setand associated spatial supports. IMPORTANT: The vertices need to be in consecutive order.
#'
#' @param pol_df a wide table format data frame identifying support of observation, or the path to a file containing the data. The data should be in wide-table format and hava column denoting the observation \code{id}, vertices \code{x1, y1, x2, y2, ...} and the time point \code{t} if required.
#' @param alpha0 sets, if needed, an averaging matrix over observations. If this specified, also a parameter \code{av_dist} needs to be specified
#' @param av_dist denotes within which distance observations are considered neighbours.
#' @param ... further arguments passed to the parent \code{Obs()} constructor for further processing.
#' @return Object of class \code{Obs_poly} (which inherits from class \code{Obs} and is thus also an observation block)
#' @keywords Observations, change of support, block
#' @export
#' @examples
#' # Create a polygon 'footprint'
#' pol_df <- data.frame(id=1,x1=0,x2=0,x3=1,x4=1,y1=0,y2=1,y3=1,y4=0,t=0)
#' pol_df <- rbind(pol_df,data.frame(id=2,x1=-0.5,x2=-0.5,x3=0.5,x4=0.5,y1=-0.5,y2=0.5,y3=0.5,y4=-0.5,t=0))
#' df <- data.frame(id=1,x=0.5,y=0.5,z=1,std=1,t=0)
#' df <- rbind(df,data.frame(id=2,x=0,y=0,z=0.2,std=0.2,t=0))
#' O <- Obs_poly(df=df,pol_df=pol_df)
#' plot(O,"z")
Obs_poly <- function(pol_df,name="Obs_poly",alpha0=NA,av_dist=NA,...)  {
   return( new("Obs_poly",name=name,pol_df,alpha0=alpha0,av_dist=av_dist,...))
}


#' @title List of links
#'
#' @description This function initialises an object of class \code{link_list} which is simply a list in which each object is restricted to be of class \code{link}.
#' @param l a list of objects of class \code{link}.
#' @return Object of class \code{link_list} (which inherits from class \code{list})
#' @keywords Observations, incidence matrix
#' @export
#' @examples
#' \dontrun{
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#'
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' L1 <- link(SURF,icesat_obs)
#' e <- link_list(list(L1))
#' }
link_list <- function(l=NULL) {
  if(!is.null(l))
    stopifnot(all(sapply(l,function(s) is(s,"link"))))
  return(new("link_list",l=l))
}

#' @title List of blocks
#'
#' @description This function initialises an object of class \code{block_list} which is simply a list in which each object is restricted to be of class \code{block}. A \code{block}
#' is also either of class \code{Obs} or class \code{process}.
#' @param l a list of objects of class \code{block}.
#' @return Object of class \code{block_list} (which inherits from class \code{list})
#' @keywords Process blocks, observation blocks
#' @export
#' @examples
#' \dontrun{
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#'
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' v <- block_list(list(O = icesat_obs, G = SURF))
#' }
block_list <- function(l=NULL) {
  if(!is.null(l))
    stopifnot(all(sapply(l,function(s) is(s,"block"))))
  return(new("block_list",l=l))
}

#' @title Graph
#'
#' @description A graph is a collection of links (collected using \code{link_list}) and blocks (collected using \code{block_list}).
#' @param e an object of class \code{link_list} containing all the edges in the graph.
#' @param v an object of class \code{block_list} containing all the blocks in the graph.
#' @return Object of class \code{Graph}
#' @keywords Process blocks, observation blocks, graph
#' @export
#' @examples
#' \dontrun{
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#'
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' L1 <- link(SURF,icesat_obs)
#' e <- link_list(list(L1))
#' v <- block_list(list(O = icesat_obs, G = SURF))
#' G <- new("Graph",e=e,v=v)
#' }
Graph <- function(e = new("link_list"),v = new("block_list")) {
  return(new("Graph",e=e,v=v))
}

#' @title Initialise a GRBF basis
#'
#' @description This function initialises an object of class \code{GRBFBasis} which defines a set of radial basis functions at pre-specified locations in 2-D
#'
#' @param x x-coordinate of GRBF centroid
#' @param y y-coordinate of GRBF centroid
#' @param std the 'length' (in terms of sigma) of the GRBF
#' @param nx the number of columns of the GRBF array
#' @param ny the number of rows of the GRBF array
#' @return Object of class \code{GRBFBasis}
#' @keywords GRBF, basis functions
#' @export
#' @examples
#' Basis <-  initGRBFbasis(x = c(0,1), y = c(0,1), std=0.1,nx=9,ny=1)
initGRBFbasis = function(x,y,std,nx,ny) {
  knots_x <- seq(x[1],x[2],length=(nx+2))
  knots_y <- seq(y[1],y[2],length=(ny+2))
  centres <- expand.grid(knots_x[-c(1,(nx+2))],knots_y[-c(1,(ny+2))])
  n <- nrow(centres)
  stds <- rep(std,n)

  fn <- pars <- list()
  for (i in 1:n) {
    fn[[i]] <-  function(pars,s) {
      return(GRBF(matrix(as.numeric(pars$centres),1,2),pars$stds,s))
    }
    pars[[i]] <- list(centres = as.matrix(centres[i,]), stds=stds[i])
  }
  df <- data.frame(x = centres[,1],
                   y = centres[,2],
                   n = 1:nrow(centres))
  pars$vars <- df
  this_basis <- new("GRBFBasis", pars=pars, n=nrow(centres), fn=fn)
  return(this_basis)
}

auto_basis <- function(m = plane(),data,nres=2,prune=1.0,type="Gaussian") {
    coords <- coordinates(data)
    bndary_seg = inla.nonconvex.hull(coords,convex=-0.05)
    loc <- scale <- NULL
    G <- list()

    xrange <- range(coords[,1])
    yrange <- range(coords[,2])

    #         ## Find possible grid points for basis locations
    for(i in 1:nres) {
        if(is(m,"plane")) {
            ## Generate mesh and use these as centres
            this_res_locs <- inla.mesh.2d(loc = matrix(apply(coords,2,mean),nrow=1),
                                          boundary = list(bndary_seg),
                                          max.edge = max(diff(xrange),diff(yrange))/(2*2.5^(i-1)),
                                          cutoff = max(diff(xrange),diff(yrange))/(3*2.5^(i-1)))$loc[,1:2]
        }  else if(is(m,"real_line")) {
            this_res_locs <- matrix(seq(xrange[1],xrange[2],length=i*6))
        } else if(is(m,"sphere")) {
            load(system.file("extdata","isea3h.rda", package = "FRK"))
            this_res_locs <- as.matrix(filter(isea3h,centroid==1,res==i)[c("lon","lat")])
        }
        ## Set scales: To 1.5x the distance to nearest basis
        ## Refine: Remove basis which are not influenced by data and re-find the scales
        for(j in 1:2) {
            D <- FRK::distance(m,this_res_locs,this_res_locs)
            diag(D) <- Inf
            this_res_scales <- apply(D,1,min)
            this_res_basis <- radial_basis(manifold = m,
                                           loc=this_res_locs,
                                           scale=ifelse(type=="Gaussian",1,1.5)*this_res_scales,
                                           type=type)
            if(j==1) {
                rm_idx <- which(colSums(eval_basis(this_res_basis,coords)) < prune)
                if(length(rm_idx) >0) this_res_locs <- this_res_locs[-rm_idx,,drop=FALSE]
            }
        }
        print(paste0("Number of basis at resolution ",i," = ",nrow(this_res_locs)))

        G[[i]] <-  radial_basis(manifold = m,
                                loc=this_res_locs,
                                scale=ifelse(type=="Gaussian",1,1.5)*this_res_scales,
                                type=type)
        G[[i]]@df$res=i
    }

    G_basis <- Reduce("concat",G)
    if(G_basis@n > nrow(data)) warning("More basis functions than data points")
    G_basis
}

radial_basis <- function(manifold=sphere(),loc=matrix(c(1,0),nrow=1),scale=1,type="Gaussian") {
    stopifnot(is.matrix(loc))
    stopifnot(dimensions(manifold) == ncol(loc))
    stopifnot(length(scale) == nrow(loc))
    stopifnot(type %in% c("Gaussian","bisquare"))
    n <- nrow(loc)
    colnames(loc) <- c(outer("loc",1:ncol(loc),FUN = paste0))

    fn <- pars <- list()
    for (i in 1:n) {
        if(type=="Gaussian") {
            fn[[i]] <-  GRBF_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        } else {
            fn[[i]] <-  bisquare_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        }
        pars[[i]] <- list(loc = matrix(loc[i,],nrow=1), scale=scale[i])
    }
    df <- data.frame(loc,scale,res=1)
    this_basis <- new("Basis", manifold=manifold, pars=pars, n=n, fn=fn, df=df)
    return(this_basis)
}

#' @title Initialise a finite element basis
#'
#' @description This function initialises an object of class \code{FEBasis} which defines a set of 'radial `tent' basis functions over a pre-specified triangulation in 2-D
#'
#' @param p \code{n} \eqn{\times} 2 matrix of vertex locations.
#' @param t \code{m} \eqn{\times} 3 matrix of triangulations. Each row identifies which rows of \code{p} make up each triangle.
#' @param M \code{n} \eqn{\times} \code{n} mass matrix: \eqn{\langle \phi, \phi^T \rangle}.
#' @param K \code{n} \eqn{\times} \code{n} stiffness matrix: \eqn{\langle \nabla\phi, \nabla\phi^T \rangle}.
#' @return Object of class \code{FEBasis}
#' @keywords finite elements, basis functions
#' @export
#' @examples
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
initFEbasis = function(p,t,M,K) {
  fn <- pars <- list()
  pars$p <- p
  pars$t <- t
  pars$M <- M
  pars$K <- K
  df <- data.frame(x = pars$p[,1],
                   y = pars$p[,2],
                   n = 1:nrow(p))
  pars$vars <- df
  # Do tessellation
  Voronoi <- deldir(pars$p[,1],
                   pars$p[,2],
                   plotit='F',
                   sort=F,
                   rw=c(min(pars$p[,1])-0.00001,
                        max(pars$p[,1])+0.00001,
                        min(pars$p[,2])-0.00001,
                        max(pars$p[,2])+.00001))
  pars$pol <- PolygonfromVoronoi(Voronoi,pars$p)

  pars$vars$area_tess_km2 = rep(0,nrow(p))
  for (i in 1:nrow(p)) {
    pars$vars$area_tess_km2[i] <- area.poly(pars$pol[[i]])
  }
  this_basis <- new("FEBasis", pars=pars, n=nrow(p), fn=fn)
  return(this_basis)}


setMethod("initialize",signature(.Object = "block"),  function(.Object,uid=NULL) {
  if(is.null(uid)) {
    .Object@uid <- round(runif(1)*1e20)
  } else {
    .Object@uid <- uid
  }
  return(.Object)
})

setMethod("initialize",signature(.Object = "GMRF"),
          function(.Object,
                   mu=matrix(0,2,1),
                   Q = sparseMatrix(i=c(1,2),j=c(1,2),x=4),
                   intrinsic=0,n=NULL,t_axis=0,
                   rep=data.frame(),name="none") {
            .Object@mu = mu
            .Object@Q = Q
            .Object@intrinsic = intrinsic
            .Object@t_axis = t_axis
            if(is.null(n)) {
              .Object@n= nrow(mu)
            } else {
              .Object@n= n
            }

            if(empty(rep)) rep <- data.frame(name=rep(name,nrow(Q)),t=rep(NA,nrow(Q)))
            .Object@rep = rep
            callNextMethod(.Object)   # to initialise block (uid)
          })


setMethod("initialize",signature(.Object = "GMRF_basis"),
          function(.Object,G=new("GMRF"),Basis=new("Basis")){
            .Object@G = G
            .Object@Basis = Basis
            if (nrow(G@rep) ==  nrow(Basis@pars$vars)) # pure spatial process
            {
              .Object@G@rep <- cbind(G@rep,Basis@pars$vars)
            } else {
              nvars <- length(which(G@rep$name == unique(G@rep$name)[1]))/nrow(Basis)/length(t_axis)
              varnum <- rep(expand.grid(1:nrow(Basis@pars$vars),1:nvars)[,2],length(t_axis))
              varnum <- c(varnum,rep(NA,nrow(G@rep) -nrow(G)))
              if(nrow(G@rep) %% nrow(Basis@pars$vars) > 0) {
                # NEEDS TO BE FIXED for when we have more covariates than 1 frame
                warning("Basis and GMRF not integral units of each other's dimensions (have covariates?), only merging first n frames")
                nframes <- floor(nrow(G@rep) / nrow(Basis@pars$vars))
                .Object@G@rep <- cbind(getDf(G)[1:(nframes*nrow(Basis)),],getDf(Basis),data.frame(varnum=varnum[1:nrow(G)]))
                extra_items <- nrow(getDf(G)) %% nrow(getDf(Basis))
                .Object@G@rep <- rbind.fill(getDf(.Object@G),tail(getDf(G),extra_items))
              } else {
                .Object@G@rep <- cbind(getDf(G),getDf(Basis),data.frame(varnum=varnum))
              }
            }
            callNextMethod(.Object,uid=.Object@G@uid)   # to initialise block (uid)
          })

setMethod("initialize",signature(.Object="VAR_Gauss"),
          function(.Object,
                   mu_fun = function(k) return(matrix(0,2,5)),
                   A_fun = function(k) return(Imat(2)),
                   B_fun =  function(k) return(emptySp()),
                   Qw_fun = function(k) return(Imat(2)),
                   t_axis = c(0:6),
                   Qb = emptySp(),
                   name="none") {

            .Object@A_fun <- A_fun
            .Object@B_fun <- B_fun
            .Object@Qw_fun <- Qw_fun
            .Object@t_axis <- t_axis
            .Object@Qb <- Qb
            n <- nrow(Qw_fun(0))
            Tn = length(t_axis)

            # My big AQA matrix
            Qw <- Qw_fun(0)
            A <- A_fun(0)
            Q_full <- Build_AQA(Qw,A,Tn)
            if (!empty(B_fun(0))) {  # if we have a B part
              B_for_sum <- vector("list",Tn)
              for(i in 0:(Tn-1)) {
                if(i == 0) {
                  #Q_beta_part <- t(A) %*% Qw %*% B_fun(i+1) - Qw %*% B_fun(i) # v1
                  #Q_beta_part <- t(A) %*% Qw %*% B_fun(i+1) - (Imat(n) - A %*% A) %*% Qw %*% B_fun(i) # v2
                  Q_beta_part <- t(A) %*% Qw %*% B_fun(i+1) - (Qw - A %*% Qw %*% A) %*% B_fun(i) # v2
                } else if (i == (Tn-1)) {
                  Q_beta_part <- rBind(Q_beta_part,- Qw %*% B_fun(i))
                }  else {
                  Q_beta_part <- rBind(Q_beta_part,t(A) %*% Qw%*% B_fun(i+1) - Qw %*% B_fun(i))
                }
                B_for_sum[[i+1]] <- t(B_fun(i)) %*% Qw %*% B_fun(i)
              }
              B_for_sum[[1]] <- t(B_fun(0)) %*% (Qw - A %*% Qw %*% A) %*% B_fun(0) # v2
              Q_full <- cBind(Q_full,Q_beta_part)
              Q_full <- rBind(Q_full,cBind(t(Q_beta_part),Reduce("+", B_for_sum) + Qb))
            }
            n_tot <- n*Tn
            Big_Q <- Q_full
            mu <- sapply(0:(Tn-1),mu_fun)
            Big_mu <- matrix(mu,ncol=1)

            field_names <- c(rep(name,n*Tn))
            if (nrow(Qb)>0)
              if (nrow(Qb) %% n == 0) {  # Assume each block of size n is one field
                  for(i in 1:(nrow(Qb)/n)) {
                    field_names <- c(field_names, rep(paste(name,i,sep=""),n))
                  }
              } else  {
                for(i in 1:(nrow(Qb))) { # otherwise assign a name to each covariate element
                  field_names <- c(field_names, paste(name,i,sep=""))
                }
            }
            rep <- data.frame(t = c(kronecker(t_axis,rep(1,n)),rep(NA,nrow(Qb))),
                              name = field_names)
            Big_mu <- rbind(Big_mu,matrix(rep(0,nrow(Qb))))

            # Cater for multi-variate fields
            #Big_Q <- as(kronecker(solve(cov_inter),Big_Q),"dgCMatrix")
            #Big_mu <- matrix(rep(mu,nrow(cov_inter)))
            #n_tot <- n_tot*nrow(cov_inter)

            callNextMethod(.Object,Big_mu,Big_Q,intrinsic=0,n=n_tot,rep=rep,t_axis=t_axis)
          })



## setMethod("initialize",signature="Obs",function(.Object,...) {

##   args<-list(...)
##   .Object@args <- args

##   if("path" %in% names(args)) {
##     cat(paste("Loading from",args$path),sep="\n")
##     data_df <- read.table(args$path,header=T)
##   } else {
##     data_df <- args$df
##   }
##   .Object@df <- data_df

##   .Object <- preprocess_obs(.Object,...)


##   if("name" %in% names(args)) {
##     .Object@df$obs_name <- as.factor(args$name)
##   }

##   if("remove_cross_ins" %in% names(args)) {
##     .Object@df <- subset(.Object@df,sqrt(x^2 + y^2) > args$remove_cross_ins)
##   }
##   if("pol" %in% names(args)) {
##     poly_points <- args$pol
##     if (!("id" %in% names(.Object@df))) stop("Need to merge by id field which is not supplied")
##     .Object@df <- merge(poly_points,.Object@df,by=c("id","t"))
##     .Object@df <- arrange(.Object@df,id,t)
##     .Object@df2 <- .expand_poly(.Object@df)
##   }
##   .Object@df$n <- 1:nrow(.Object@df)
##   .Object@n <- nrow(.Object@df)

##   if("cmweq" %in% names(.Object@df)) {
##     .Object@df$z <-  as.numeric(.Object@df$cmweq)*0.01*.Object@df$area2    #Convert to Mt
##     .Object@df2$z <-  as.numeric(.Object@df2$cmweq)*0.01*.Object@df2$area2    #Convert to Mt
##     .Object@df$std <-  as.numeric(.Object@df$std)*0.01*.Object@df$area2    #Convert to Mt
##     .Object@df2$std <-  as.numeric(.Object@df2$std)*0.01*.Object@df2$area2    #Convert to Mt
##   }

##   if("alpha0" %in% names(args)) {
##     if(!("av_dist" %in% names(args))) stop("Cannot specify alpha0 without averaging distance")
##     .Object@args$P <-  Find_Smooth_mat(subset(.Object@df,t==0),args$alpha0,args$av_dist)
##   }

##   callNextMethod(.Object)} )

# ... is passed on to preprocess_obs
setMethod("initialize",signature="Obs",function(.Object,df,name=NA,remove_cross_ins=0,pol=NA,alpha0=NA,av_dist=NA,...) {

    args<-list(...)
    args <- c(args,df=df,name=name,remove_cross_ins=remove_cross_ins,pol=pol,alpha0=alpha0,av_dist=av_dist)
    .Object@args <- args

    stopifnot((is.character(df)) | is.data.frame(df))

    if(is.character(df)) {
        cat(paste("Loading from",df),sep="\n")
        data_df <- read.table(df,header=T)
    } else {
        data_df <- df
    }

    .Object@df <- data_df
    .Object <- preprocess_obs(.Object,...)
    if (is.null(data_df$obs_name))
      .Object["obs_name"] <- as.factor(name)
    if(remove_cross_ins > 0) {
        .Object@df <- subset(.Object@df,sqrt(x^2 + y^2) > remove_cross_ins)
    }

    if(!is.na(pol[1])) {
        poly_points <- pol
        if (!("id" %in% names(.Object@df))) stop("Need to merge by id field which is not supplied")
        .Object@df <- merge(poly_points,.Object@df,by=c("id","t"))
        .Object@df <- arrange(.Object@df,id,t)
        .Object@df2 <- .expand_poly(.Object@df)
    }
    .Object@df$n <- 1:nrow(.Object@df)
    .Object@n <- nrow(.Object@df)

    if("cmweq" %in% names(.Object@df)) {
        .Object@df$z <-  as.numeric(.Object@df$cmweq)*0.01*.Object@df$area2    #Convert to Mt
        .Object@df2$z <-  as.numeric(.Object@df2$cmweq)*0.01*.Object@df2$area2    #Convert to Mt
        .Object@df$std <-  as.numeric(.Object@df$std)*0.01*.Object@df$area2    #Convert to Mt
        .Object@df2$std <-  as.numeric(.Object@df2$std)*0.01*.Object@df2$area2    #Convert to Mt
    }


    if(!is.na(alpha0)) {
        if(is.na(alpha0)) stop("Cannot specify alpha0 without averaging distance")
        .Object@args$P <-  Find_Smooth_mat(subset(.Object@df,t==0),alpha0,av_dist)
    }

    callNextMethod(.Object)})



## setMethod("initialize",signature="Obs_poly",function(.Object,...) {
##   args <- list(...)
##   cat("Forming polygons from indices. Please wait...",sep="\n")
##   if ("poly_path" %in% names(args)) {
##     pol <- read.table(args$poly_path,header=T)
##   } else if ("pol_df" %in% names(args))   {
##     pol <- args$pol_df
##   } else {
##     stop("No arguments")
##   }
##   x_ind <- grep("x[[:digit:]]",names(pol))
##   y_ind <- grep("y[[:digit:]]",names(pol))
##   x <- as.matrix(pol[,x_ind]);
##   y <- as.matrix(pol[,y_ind]);
##   pol_list <- pol_from_xy(x,y,order_vertices=F)
##   .Object@pol <- pol_list
##   .Object@df2 <- data.frame()
##   callNextMethod(.Object,pol=pol,...)
## })

setMethod("initialize",signature="Obs_poly",function(.Object,pol_df,...) {
  cat("Forming polygons from indices. Please wait...",sep="\n")
  stopifnot((is.character(pol_df)) | is.data.frame(pol_df))
  if(is.character(pol_df)) {
    pol <- read.table(pol_df,header=T)
  } else {
    pol <- pol_df
  }
  x_ind <- grep("x[[:digit:]]",names(pol))
  y_ind <- grep("y[[:digit:]]",names(pol))
  x <- as.matrix(pol[,x_ind]);
  y <- as.matrix(pol[,y_ind]);
  pol_list <- pol_from_xy(x,y,order_vertices=F)
  .Object@pol <- pol_list
  .Object@df2 <- data.frame()
  callNextMethod(.Object,pol=pol,...)
})



setMethod("initialize",signature(.Object="link"), function(.Object,from=new("block"),to=new("block")) {
  .Object@from=from
  .Object@to=to
  return(.Object)
})

setMethod("initialize",signature(.Object = "linkGO"),  function(.Object,from=new("process"),to=new("Obs"),
                                                                n_grid = NULL,Cmat = NULL, mul_factor = NULL,
                                                                mulfun = NULL, mask = NULL, muldata = NULL,
                                                                md5_wrapper = NULL) {

  if(!is.null(Cmat)) {
      .Object@Cmat <- Cmat
  } else {
    if(class(from@Basis) == "FEBasis") {
      t_axis <- from@G@t_axis
      Tn <- length(t_axis)
      n <- from@Basis@n # number of elements
      if(class(from@G) == "VAR_Gauss") {
        #Tn <- from@G@Tn # if spatio-temporal, Tn is the number of time-points we need to infer at
        nvariates <- (from@G@n/Tn) / n
      } else {
        #Tn <- max(to@df$t)+1 # if spatial, Tn is the number of time points we have observations of...
        nvariates <- 1
      }

      if(!(round(nvariates) == nvariates)) stop("Stopping: Cannot calculate number of variates. Are they on the same basis?")
      if(nvariates < 1) stop("Stopping: Cannot calculate number of variates")
      if(nvariates > 1) cat("Multi-variate system detected",sep="\n")
      if (!all(is.na(getDf(from)$t))) # If not a spatial process
        if(!all(unique(getDf(to)$t) %in% unique(getDf(from)$t))) stop("Cannot have observations which exceed the modelling boundary in t_axis")
      if(!(is.null(mul_factor))) {
        if(!(length(mul_factor) == nvariates)) stop("Need as many mul_factors as variates")
      } else {
        mul_factor <- rep(1,nvariates)
      }

      Cmats <- vector("list",Tn)

      if(any(diff(to@df$t) < 0)) stop("Can only do link using data which is ordered temporally. Please order the data and redo.")

      stopifnot(all(diff(to@df$t) >= 0))
      for(i in seq_along(t_axis)) {
        to_sub <- to
        to_sub@df <- subset(to@df,t==t_axis[i]) # find data points at this time point
        if(class(to_sub) == "Obs_poly") {
          to_sub@pol <- to_sub@pol[which(to_sub@df$t==t_axis[i])] # if polygon also extract polygons
        }
        if(nrow(to_sub)==0) {  # if no data points at this time point
          Cmats[[i]] <- matrix(0,0,n)  # create empty matrix
        } else {
          Cmats[[i]] <- .find_inc_matrix(from@Basis,to_sub,mulfun = mulfun, muldata = muldata, mask = mask, n_grid = n_grid,
                                         md5_wrapper = md5_wrapper) # otherwise compute the matrix
        }
        C <- Cmats[[i]]  # store incidence matrix at this time point

        j <- 1
        Cmats[[i]] <- matrix(0,nrow(C),0) # if we have more than one variate, concatenate another version of inc. matrix
        while (j <= nvariates) {
          if(mul_factor[j] == 0)  {
            C2 <- Zeromat(nrow(C),ncol(C))
          } else {
            C2 <- mul_factor[j] * C
          }
          Cmats[[i]]  <- cBind(Cmats[[i]],C2)
          j <- j+1
        }
      }



      if(class(from@G) == "VAR_Gauss") { # if spatio-temporal process
        .Object@Cmat <- do.call("bdiag",Cmats)  # diagonally bind matricies
        .Object@Cmat <- cBind(.Object@Cmat,Zeromat(nrow(.Object@Cmat),nrow(from@G@Qb))) # and add on covariate effect
      } else { # if spatial process
        .Object@Cmat <- do.call("rBind",Cmats) # vertically bind matrices
        # PS: No Qb because this is a repeatedly observed spatial field

      }

    }
  }
  callNextMethod(.Object,from,to)
})


setMethod("initialize",signature(.Object = "linkGG"),  function(.Object,from=new("process"),to=new("process"),cov_inter = matrix(0,0,0)) {

  #stop("Links not specified between processes. Use the cov_inter option when definint processes for coregionalisation.")
  .Object@cov_inter = cov_inter
  callNextMethod(.Object,from,to)
})

setMethod("initialize",signature(.Object="link_list"),function(.Object,l=NULL){
  if(is.null(l)) {
    .Object@.Data = list(L1 = new("link"))
  } else {
    .Object@.Data = l
  }
  return(.Object)})
setMethod("initialize",signature(.Object="block_list"),  function(.Object,l=NULL){
  if(is.null(l)) {
    .Object@.Data = list(G=new("GMRF"),O=new("Obs"))
  } else {
    .Object@.Data = l}
  return(.Object)})

