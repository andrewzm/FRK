#' @title Construct a set of local basis functions
#' @description Construct a set of local basis functions based on pre-specified location and scale parameters.
#' @param manifold object of class \code{manifold}, for example, \code{sphere}
#' @param loc a matrix of size \code{n} by \code{dimensions(manifold)} indicating centres of basis functions
#' @param scale vector of length \code{n} containing the scale parameters of the basis functions; see details
#' @param type either ``Gaussian'', ``bisquare,'' ``exp,'' or ``Matern32''
#' @details This functions lays out local basis functions in a domain of interest based on pre-specified location and scale parameters. If \code{type} is ``Gaussian'', then
#' \deqn{\phi(u) = \exp\left(-\frac{\|u \|^2}{2\sigma^2}\right),}
#' and \code{scale} is given by \eqn{\sigma}, the standard deviation. If \code{type} is ``bisquare'', then
#'\deqn{\phi(u) = \left(1- \left(\frac{\| u \|}{R}\right)^2\right)^2 I(\|u\| < R),}
#' and \code{scale} is given by \eqn{R}, the range of support of the bisquare function. If the \code{type} is ``exp'', then
#'\deqn{\phi(u) = \exp\left(-\frac{\|u\|}{ \tau}\right),}
#' and \code{scale} is given by \eqn{\tau}, the e-folding length. If \code{type} is ``Matern32'', then
#'\deqn{\phi(u) = \left(1 + \frac{\sqrt{3}\|u\|}{\kappa}\right)\exp\left(-\frac{\sqrt{3}\| u \|}{\kappa}\right),}
#' and \code{scale} is given by \eqn{\kappa}, the function's scale.
#' @examples
#' library(ggplot2)
#' G <-  local_basis(manifold = real_line(),
#'                    loc=matrix(1:10,10,1),
#'                    scale=rep(2,10),
#'                    type="bisquare")
#' show_basis(G)
#' @export
local_basis <- function(manifold=sphere(),loc=matrix(c(1,0),nrow=1),scale=1,type="Gaussian") {
    stopifnot(is.matrix(loc))
    stopifnot(dimensions(manifold) == ncol(loc))
    stopifnot(length(scale) == nrow(loc))
    stopifnot(type %in% c("Gaussian","bisquare","exp","Matern32"))
    n <- nrow(loc)
    colnames(loc) <- c(outer("loc",1:ncol(loc),FUN = paste0))

    fn <- pars <- list()
    for (i in 1:n) {
        if(type=="Gaussian") {
            fn[[i]] <-  .GRBF_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        } else if (type=="bisquare") {
            fn[[i]] <-  .bisquare_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        } else if (type=="exp") {
            fn[[i]] <-  .exp_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        } else if (type=="Matern32") {
            fn[[i]] <-  .Matern32_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        }

        pars[[i]] <- list(loc = matrix(loc[i,],nrow=1), scale=scale[i])
    }
    df <- data.frame(loc,scale,res=1)
    this_basis <- new("Basis", manifold=manifold, pars=pars, n=n, fn=fn, df=df)
    return(this_basis)
}

#' @title Automatic basis-function placement
#' @description Generate automatically a set of local basis functions in the domain, and automatically prune in regions of sparse data.
#' @param manifold object of class \code{manifold}, for example, \code{sphere} or \code{plane}
#' @param data object of class \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame} containing the data on which basis-function placement is based, or a list of these; see details
#' @param regular an integer indicating the number of regularly-placed basis functions at the first resolution. In two dimensions, this dictates smallest number of basis functions in a row or column at the lowest resolution. If \code{regular=0}, an irregular grid is used, one that is based on the triangulation of the domain with increased mesh density in areas of high data density, see details
#' @param nres if \code{manifold = real_line()} or \code{manifold = plane()}, then \code{nres} is the number of basis-function resolutions to use. If \code{manifold = sphere()}, then \code{nres} is the resolution number of the ISEA3H grid to use and and can also be a vector indicating multiple resolutions
#' @param prune a threshold parameter which dictates when a basis function is considered irrelevent or unidentifiable, and thus removed, see details
#' @param max_basis maximum number of basis functions. This overrides the parameter \code{nres}
#' @param subsamp the maximum amount of data points to consider when carrying out basis-function placement: these data objects are randomly sampled from the full dataset. Keep this number fairly high (on the order of 10^5) otherwise high resolution basis functions may be spuriously removed
#' @param type the type of basis functions to use; see details
#' @param isea3h_lo if \code{manifold = sphere()}, this argument dictates which ISEA3H resolution is the lowest one that should be used for basis-function placement
#' @param bndary a \code{matrix} containing points containing the boundary. If \code{regular == 0} this can be used to define a boundary in which irregularly-spaced basis functions are placed
#' @param verbose a logical variable indicating whether to output a summary of the basis functions created or not
#' @param ... unused
#' @details This function automatically places basis functions within the domain of interest. If the domain is a plane or the real line, then the object \code{data} is used to establish the domain boundary.
#'
#'
#'The argument \code{type} can be either ``Gaussian'' in which case
#'\deqn{\phi(u) = \exp\left(-\frac{\|u \|^2}{2\sigma^2}\right),}
#'``bisquare'', in which case
#'\deqn{\phi(u) = \left(1- \left(\frac{\| u \|}{R}\right)^2\right)^2 I(\|u\| < R),}
#'``exp'', in which case
#'\deqn{\phi(u) = \exp\left(-\frac{\|u\|}{ \tau}\right),}
#' or ``Matern32'', in which case
#'\deqn{\phi(u) = \left(1 + \frac{\sqrt{3}\|u\|}{\kappa}\right)\exp\left(-\frac{\sqrt{3}\| u \|}{\kappa}\right),}
#' The parameters \eqn{\sigma, R, \tau} and \eqn{\kappa} are \code{scale} arguments.
#'
#' If the manifold is the real line, the basis functions are placed regularly inside the domain, and the number of basis functions at the lowest resolution is dictated by the integer parameter \code{regular} which has to be greater than zero. On the real line, each subsequent resolution has twice as many basis functions. The scale of the basis function is set based on the minimum distance between the centre locations following placement. The scale is equal to the minimum distance if the type of basis function is Gaussian, exponential or Matern32, and is equal to 1.5 times this value if the function is bisquare.
#'
#' If the manifold is a plane, and \code{regular > 0}, then basis functions are placed regularly within the bounding box of \code{data}, with the smallest number of basis functions in each row or column equal to the value of \code{regular} in the lowest resolution. Subsequent resolutions have twice the number of basis functions in each row or column. If \code{regular = 0}, then the function \code{INLA::inla.nonconvex.hull} is used to construct a (non-convex) hull around the data. The buffer and smoothness of the hull is determined by the parameter \code{convex}. Once the domain boundary is found,  \code{INLA::inla.mesh.2d} is used to construct a triangular mesh such that the node vertices coincide with data locations, subject to some minimum and maximum triangular side length constraints. The result is a mesh which is dense in regions of high data density and not dense in regions of sparse data. Even in this case, the scale is taken to be a function of the minimum distance between basis function centres, as detailed above. This may be changed in a future revision.
#'
#' If the manifold is the surface of a sphere, then basis functions are placed on the centroids of the discrete global grid (DGG), with the first basis resolution corresponding to the first resolution of the DGG (ISEA3H resolution 0, which yields 12 basis functions globally).  It is not recommended to go above \code{nres == 5} (ISEA3H resolutions 0--4) for the whole sphere, which would yield a total of 1220 basis functions. Up to ISEA3H resolution 6 is available with \code{FRK}; for higher resolutions please install \code{dggrids} from \code{https://github.com/andrewzm/dggrids}.
#'
#' Basis functions that are not influenced by data points may hinder convergence of the EM algorithm, since the associated hidden states are by and large unidentifiable. We hence provide a means to automatically remove such basis functions through the parameter \code{prune}. The final set only contains basis functions for which the column sums in the associated matrix \eqn{S} (which, recall, is the value/average of the basis functions at/over the data points/polygons) is greater than \code{prune}. If \code{prune == 0}, no basis functions are removed from the original design.
#' @examples
#' library(sp)
#' library(ggplot2)
#'
#' ### Create a synthetic dataset
#' d <- data.frame(lon = runif(n=1000,min = -179, max = 179),
#'                 lat = runif(n=1000,min = -90, max = 90),
#'                 z = rnorm(5000))
#' coordinates(d) <- ~lon + lat
#' proj4string(d)=CRS("+proj=longlat +ellps=sphere")
#'
#' ### Now create basis functions over sphere
#' G <- auto_basis(manifold = sphere(),data=d,
#'                 nres = 2,prune=15,
#'                 type = "bisquare",
#'                 subsamp = 20000)
#'
#' ### Plot
#' show_basis(G,draw_world())
#'
#' @export
auto_basis <- function(manifold = plane(),
                       data,
                       regular=1,
                       nres=2,
                       prune=0,
                       max_basis = NULL,
                       subsamp=10000,
                       type="bisquare",
                       isea3h_lo = 0,
                       bndary = NULL,
                       verbose = 0L,
                       ...) {
    m <- manifold
    if(!is(m,"manifold"))
        stop("manifold needs to be an object of class manifold")
    if(!is.numeric(prune) | prune < 0)
        stop("prune needs to be greater than zero")
    if(!is.numeric(subsamp) | subsamp < 0)
        stop("subsamp needs to be greater than zero")
    if(!type %in% c("bisquare","Gaussian","exp","Matern32"))
        stop("type of basis functions must be 'Gaussian,' 'bisquare,' 'exp,' or 'Matern32'.")
    if((is(m,"sphere")  | is(m,"real_line")) & regular == 0)
        stop("Irregular basis only available on planes")
    if(!(is(isea3h_lo,"numeric")))
        stop("isea3h_lo needs to be an integer greater than 0")
    if(!(is.numeric(nres) | is.null(nres)))
        stop("nres needs to be greater than zero or NULL")
    if(!is.null(max_basis)) {
        print("...Automatically choosing functions...")
        tot_basis <- 0
        tot_data <- length(data)
        nres <- 1
        while(tot_basis <= max_basis) {
            nres <- nres + 1
            G <- .auto_basis(manifold =manifold,
                            data=data,
                            prune =0,regular=regular,nres=nres,
                            subsamp=subsamp,type=type,isea3h_lo = isea3h_lo,
                            bndary=bndary, verbose=0)
            tot_basis <- nbasis(G)
        }
        S <- eval_basis(G,data)
        prune <- (colSums(S)[rev(order(colSums(S)))])[round(max_basis)] + 1e-10
    }

    .auto_basis(manifold=manifold,data=data,regular=regular,nres=nres,
                prune=prune,subsamp=subsamp,type=type,isea3h_lo = isea3h_lo,
                bndary=bndary, verbose=verbose)

}


.auto_basis <- function(manifold = plane(),
                       data,
                       regular=1,
                       nres=2,
                       prune=0,
                       subsamp=10000,
                       type="Gaussian",
                       isea3h_lo = 0,
                       bndary = NULL,
                       verbose = 0L) {

    m <- manifold
    isea3h <- centroid <- res <- NULL #(suppress warnings, these are loaded from data)
    coords <- coordinates(data)

    if(is(m,"plane") & regular == 0 & is.null(bndary)) {
         if(!requireNamespace("INLA"))
             stop("For irregularly-placed basis-function generation INLA needs to be installed
                  for constructing basis function centres. Please install it
                  using install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/stable\")")
    }

    if(nrow(coords)>subsamp) {
        coords <- coords[sample(nrow(coords),size=subsamp,replace=FALSE),]
    }

    xrange <- range(coords[,1])
    yrange <- range(coords[,2])


    if(is(m,"plane") & regular == 0) {
        if(is.null(bndary)) {
            bndary_seg = INLA::inla.nonconvex.hull(coords,concave = 0)
        } else {
            if(!is(bndary,"matrix"))
                stop("bndary needs to be a matrix of points")
            bndary_seg <- INLA::inla.mesh.segment(bndary)
        }
    }

    if(is(m,"plane") & regular > 0) {
        asp_ratio <- diff(yrange) / diff(xrange)
        if(asp_ratio < 1) {
            ny <- regular
            nx <- round(ny / asp_ratio)
        } else {
            nx <- regular
            ny <- round(nx * asp_ratio)
        }
    }

    loc <- scale <- NULL
    G <- list()

    if(is(m,"sphere")) {
        isea3h <- load_dggrids(res = nres) %>%
                  dplyr::filter(res >= isea3h_lo)
        #nres <- nres - isea3h_lo
    }


    ## Find possible grid points for basis locations
    for(i in 1:nres) {
        if(is(m,"plane") & (regular == 0)) {
            ## Generate mesh and use these as centres
            this_res_locs <- INLA::inla.mesh.2d(loc = matrix(apply(coords,2,mean),nrow=1),
                                          boundary = list(bndary_seg),
                                          max.edge = max(diff(xrange),diff(yrange))/(2*2.5^(i-1)),
                                          cutoff = max(diff(xrange),diff(yrange))/(3*2.5^(i-1)))$loc[,1:2]
        } else if(is(m,"plane") & (regular> 0)) {
           xgrid <- seq(xrange[1], xrange[2], length =nx*(3^(i)))
           ygrid <- seq(yrange[1], yrange[2], length =ny*(3^(i)))
            ## Generate mesh and use these as centres
            #xgrid <- seq(xrange[1] + diff(xrange)/(2*nx*i), xrange[2] - diff(xrange)/(2*nx*i), length =nx*(3^(i-1)))
            #ygrid <- seq(yrange[1] + diff(yrange)/(2*ny*i), yrange[2] - diff(yrange)/(2*ny*i), length =ny*(3^(i-1)))
            this_res_locs <- xgrid %>%
                expand.grid(ygrid) %>%
                as.matrix()
        } else if(is(m,"real_line")) {
            this_res_locs <- matrix(seq(xrange[1],xrange[2],length=i*regular))
        } else if(is(m,"sphere")) {
            this_res_locs <- as.matrix(filter(isea3h,centroid==1,res==(i-1 + isea3h_lo))[c("lon","lat")])
        }
        ## Set scales: To 1.5x the distance to nearest basis if bisquare
        ## Refine: Remove basis which are not influenced by data and re-find the scales
        for(j in 1:2) {
            D <- FRK::distance(m,this_res_locs,this_res_locs)
            if(nrow(D) == 1) {
                this_res_scales <-max(diff(xrange),diff(yrange))/2
            } else {
                diag(D) <- Inf
                this_res_scales <- apply(D,1,min)
            }

            if(nrow(D) >0)
                this_res_basis <- local_basis(manifold = m,
                                           loc=this_res_locs,
                                           scale=ifelse(type=="bisquare",1.5,1.5)*this_res_scales,
                                           type=type)
            if(prune > 0 & nrow(D)>0) {
                if(j==1) {
                    rm_idx <- which(colSums(eval_basis(this_res_basis,coords)) < prune)
                    if(length(rm_idx) == length(this_res_scales))
                        warning("prune is too large -- all functions at a resolution removed.
                             Consider also removing number of resolutions.")
                    if(length(rm_idx) >0) this_res_locs <- this_res_locs[-rm_idx,,drop=FALSE]
                }
            } else {
                break
            }
        }
        if(verbose) print(paste0("Number of basis at resolution ",i," = ",nrow(this_res_locs)))

        if(nrow(D) > 0) {
            G[[i]] <-  local_basis(manifold = m,
                                    loc=this_res_locs,
                                    scale=ifelse(type=="bisquare",1.5,1.5)*this_res_scales,
                                    type=type)
            G[[i]]@df$res=i
        }
    }

    G_basis <- Reduce("concat",G)
    if(G_basis@n > nrow(data)) warning("More basis functions than data points")
    if(G_basis@n > 2000) warning("More than 2000 basis functions")
    G_basis
}

#' @title Add the time coordinate to 2D spatial basis functions
#' @description Given a set of 2D spatial basis functions and a vector of knots in time, this function repeats the spatial basis at every temporal knot, adding the third dimension (i.e., time) to the centroid as appropriate.
#' @param G_spatial an object of class Basis on a 2D manifold
#' @param t_knots a vector of numbers locating the knots in time
#' @param manifold a 3D space-time manifold, typically \code{STsphere()} or \code{STplane()}
#' @examples
#' G_spatial <-  local_basis(manifold = sphere(),
#'                    loc=matrix(runif(20,min=-90,max=90),10,2),
#'                    scale=rep(20,10),
#'                    type="bisquare")
#' G_space_time <- sp_to_ST_basis(G_spatial,1:10,manifold=STsphere())
#' library(ggplot2)
#' show_basis(G_space_time)
#' @export
sp_to_ST_basis <- function(G_spatial,t_knots = 1,manifold=STsphere()) {
    stopifnot(dimensions(manifold(G_spatial))==2)
    stopifnot(dimensions(manifold)==3)
    stopifnot(is.numeric(t_knots))

    n <- G_spatial@n
    G <- list()
    for(i in seq_along(t_knots)) {
        Gt <- G_spatial
        sapply(1:n, function(j) {
            this_c <-   get("c",environment(Gt@fn[[j]]))    # retrieve centroid
            new_c <- cbind(this_c,t_knots[i])               # add time coordinate
            assign("c",new_c,environment(Gt@fn[[j]]))
        })
        Gt@df <- cbind(Gt@df,loc3=t_knots[i])
        Gt@manifold <- manifold

        G[[i]] <- Gt
    }
    G <- Reduce("concat",G)
}

#' @rdname TensorP
#' @aliases TensorP,Basis-Basis-method
setMethod("TensorP",signature(Basis1="Basis",Basis2="Basis"),function(Basis1,Basis2) {
    if(nres(Basis2) > 1) stop("Only one basis can be multiresolution (Basis 1)")
    df1 <- Basis1@df
    df2 <- Basis2@df
    n1 <- dimensions(Basis1)
    n2 <- dimensions(Basis2)
    expand.grid(df1[,1:n1,drop=FALSE],df2[,1:n2,drop=FALSE])
    df <- cbind(df1[rep(1:nrow(df1),times = nrow(df2)),1:n1], ## One resolution in df2 being assumed
                df2[rep(1:nrow(df2),each = nrow(df1)),1:n2],
                df1[rep(1:nrow(df1),times = nrow(df2)),"res"])
    names(df) <- c(paste0("loc",1:(n1 + n2)),"res")

    new("TensorP_Basis",
        Basis1=Basis1,
        Basis2=Basis2,
        n = Basis1@n * Basis2@n,
        df = df)
})


#' @rdname eval_basis
#' @aliases eval_basis,Basis-matrix-method
setMethod("eval_basis",signature(basis="Basis",s="matrix"),function(basis,s,output = "matrix"){
    stopifnot(output %in% c("list","matrix"))
    space_dim <- dimensions(manifold(basis))
    if(opts_FRK$get("parallel") > 1L) {
        n <- nrow(s)
        batching=cut(1:n,breaks = seq(0,n+10000,
                                      by=10000),labels=F)

        clusterExport(opts_FRK$get("cl"),
                      c("batching","basis","s","space_dim","output"),envir=environment())
        pnt_eval_list <- parLapply(opts_FRK$get("cl"),1:max(unique(batching)),
                                  function(i) {
                                      idx <- which(batching == i)
                                      return(.point_eval_fn(basis@fn,
                                                     s[idx,1:space_dim,drop=F],output))
                                  })
        clusterEvalQ(opts_FRK$get("cl"), {gc()})
        do.call(rBind,pnt_eval_list)
    } else  {
        .point_eval_fn(basis@fn,s[,1:space_dim,drop=F],output)
    }
})

#' @rdname eval_basis
#' @aliases eval_basis,Basis-SpatialPointsDataFrame-method
setMethod("eval_basis",signature(basis="Basis",s="SpatialPointsDataFrame"),function(basis,s,output = "matrix"){
    stopifnot(output %in% c("matrix","list"))
    space_dim <- dimensions(manifold(basis))

    if(opts_FRK$get("parallel") > 1L) {
        n <- length(s)
        batching=cut(1:n,breaks = seq(0,n+10000,
                                      by=10000),labels=F)


        clusterExport(opts_FRK$get("cl"),
                      c("batching","basis","s","space_dim","output"),envir=environment())
        pnt_eval_list <- parLapply(opts_FRK$get("cl"),1:max(unique(batching)),
                                   function(i) {
                                       idx <- which(batching == i)
                                       return( .point_eval_fn(basis@fn,
                                                              coordinates(s)[idx,1:space_dim,drop=F],output))
                                   })
        clusterEvalQ(opts_FRK$get("cl"), {gc()})
        do.call(rBind,pnt_eval_list)
    } else  {
        .point_eval_fn(basis@fn,coordinates(s)[,1:space_dim,drop=F],output)
    }
})

#' @rdname eval_basis
#' @aliases eval_basis,Basis-SpatialPolygonsDataFrame-method
setMethod("eval_basis",signature(basis="Basis",s="SpatialPolygonsDataFrame"),function(basis,s,output = "matrix"){
    stopifnot(output %in% c("matrix","list"))
    X <- list()
    print("Averaging over polygons")

    if(opts_FRK$get("parallel") > 1L) {

        ## parLapply version not tested yet
        clusterExport(opts_FRK$get("cl"),
                      c("basis","s"),envir=environment())
        X <- parLapply(opts_FRK$get("cl"),1:length(s), function(i) {
            samps <- .samps_in_polygon(basis,s,i)
            colSums(.point_eval_fn(basis@fn,samps))/nrow(samps)
        })
        clusterEvalQ(opts_FRK$get("cl"), {gc()})

        # X <- mclapply(1:length(s), function(i) {
        #     samps <- .samps_in_polygon(basis,s,i)
        #     colSums(.point_eval_fn(basis@fn,samps))/nrow(samps)
        # },mc.cores = opts_FRK$get("parallel"))

    } else  {
        X <- lapply(1:length(s), function(i) {
            samps <- .samps_in_polygon(basis,s,i)
            colSums(.point_eval_fn(basis@fn,samps))/nrow(samps)
        })
    }
    X <- Reduce("rBind",X)
    as(X,"Matrix")
})

#' @rdname eval_basis
#' @aliases eval_basis,Basis-STIDF-method
setMethod("eval_basis",signature(basis="Basis",s="STIDF"),function(basis,s,output = "matrix"){
    stopifnot(output %in% c("matrix","list"))
    space_dim <- dimensions(manifold(basis))
    .point_eval_fn(basis@fn,cbind(coordinates(s),t=s@data$t)[,1:space_dim,drop=F],output)
})

#' @rdname eval_basis
#' @aliases eval_basis,TensorP_Basis-matrix-method
setMethod("eval_basis",signature(basis="TensorP_Basis",s="matrix"),function(basis,s,output = "matrix"){
    n1 <- dimensions(manifold(basis@Basis1))
    S1 <- eval_basis(basis@Basis1,s[,1:n1,drop=FALSE],output)
    S2 <- eval_basis(basis@Basis2,s[,-(1:n1),drop=FALSE],output)

#     warning("Changed to matrix, improve")
#     S <- matrix(0,nrow(S1),ncol(S1)*ncol(S2))
#     for(i in 1:ncol(S1)) {
#         XX <-  ((S1[,i] * S2) %>% as.matrix())
#         S[,((i-1)*ncol(S2)+1):(i*ncol(S2))] <- XX
#     }
    #XX <- lapply(1:ncol(S1),function(i)  (S1[,i] * S2))
    ## Order: First space then time
    XX <- lapply(1:ncol(S2),function(i)  (S2[,i] * S1))
    S <- quickcBind(XX)
    S
})


#' @rdname eval_basis
#' @aliases eval_basis,TensorP_Basis-STIDF-method
setMethod("eval_basis",signature(basis="TensorP_Basis",s = "STIDF"),function(basis,s,output = "matrix"){
    n1 <- dimensions(manifold(basis@Basis1))
    slocs <- coordinates(s)
    tlocs <- matrix(s@data$t)

    S1 <- eval_basis(basis@Basis1,slocs[,,drop=FALSE],output)
    S2 <- eval_basis(basis@Basis2,tlocs[,,drop=FALSE],output)

    XX <- lapply(1:ncol(S2),function(i)  (S2[,i] * S1))
    S <- quickcBind(XX)
    S
})


#' @rdname eval_basis
#' @aliases eval_basis,TensorP_Basis-STFDF-method
setMethod("eval_basis",signature(basis="TensorP_Basis",s = "STFDF"),function(basis,s,output = "matrix"){
    n1 <- dimensions(manifold(basis@Basis1))
    slocs <- coordinates(s)
    tlocs <- matrix(s@data$t)
    nt <- length(s@time)

    S1 <- eval_basis(basis@Basis1,s[,1],output)
    S1 <- do.call("rBind",lapply(1:nt,function(x) S1))
    S2 <- eval_basis(basis@Basis2,tlocs[,,drop=FALSE],output)
    XX <- lapply(1:ncol(S2),function(i)  (S2[,i] * S1))
    S <- quickcBind(XX)
    S
})


#' @rdname local_basis
#' @export
radial_basis <- function(manifold=sphere(),loc=matrix(c(1,0),nrow=1),scale=1,type="Gaussian") {
    stop("radial_basis is deprecated. Please use local_basis instead")

}

#' @rdname nres
#' @aliases nres_basis,Basis-method
setMethod("nres",signature(b="Basis"),function(b){ length(unique(b@df$res))})

#' @rdname nres
#' @aliases nres_basis,Basis-method
setMethod("nres",signature(b="TensorP_Basis"),function(b){nres(b@Basis1) * nres(b@Basis2)})


#' @rdname nres
#' @aliases nres_SRE,SRE-method
setMethod("nres",signature(b="SRE"),function(b){ nres(b@basis)})

setMethod("BuildD",signature(G="Basis"),function(G){
    res <- NULL # suppress bindings (it's in the data frame)
    nres <- nres(G)
    m <- manifold(G)
    D_basis = lapply(1:nres,function(i)  {
        x1 <- filter(G@df,res == i)[,1:dimensions(m)] %>% as.matrix()
        distance(m,x1,x1)
})})

setMethod("BuildD",signature(G="TensorP_Basis"),function(G){
    nres1 <- nres(G@Basis1)
    nres2 <- nres(G@Basis2)
    stopifnot(nres2 == 1)
    D_basis <- list(Basis1 = BuildD(G@Basis1),
                    Basis2 = BuildD(G@Basis2))
    D_basis
    })

.point_eval_fn <- function(flist,s,output="matrix") {

    x <- do.call("cbind",sapply(flist,function(f) f(s),simplify=FALSE))
    as(x,"Matrix")

    ## Rhipe Hadoop version (currently disabled)
    ## The below works but likely to be slower.. whole prediction should be parallelised

    # envlist <- lapply(flist, function(f) environment(f))
    # flist <- lapply(flist, function(f) parse(text = deparse(f)))
    # x <- rhwrapper(Ntot = nrow(s),
    #                N = 4000,
    #                type="Matrix",
    #                f_expr = .rhpoint_eval_fn,
    #                flist=flist,
    #                envlist = envlist,
    #                s=s)
    # as(data.matrix(x),"Matrix")

}

.samps_in_polygon <- function(basis,s,i) {
    nMC <- 1000
    if(is(basis@manifold,"plane")) {
        samps <- coordinates(spsample(s@polygons[[i]],n=nMC,type="random"))
    } else if(is(basis@manifold,"sphere")){

        ## Find coordinates
        coords <- data.frame(slot(s@polygons[[i]]@Polygons[[1]],"coords"))

        ## Find lat/lon range
        rangelon <- range(coords$lon) + 180
        rangelat <- range(coords$lat) + 90

        ## Find limits in transformed space
        rangeu <- rangelon/360
        rangev <- (cos(rangelat*2*pi/360) + 1)/2

        ## Sample in transformed space
        samps <- cbind(runif(nMC,min=min(rangeu), max=max(rangeu)),
                       runif(nMC,min=min(rangev), max=max(rangev)))

        ## Re-transform back to spherical coordinates
        samps[,1] <- 360*samps[,1] - 180
        samps[,2] <- acos(2*samps[,2] -1) * 360 / (2*pi) - 90

        ## Find which points are in polygon (approx. 70%)
        pip <- over(SpatialPoints(samps),
                    SpatialPolygons(list(s@polygons[[i]]),1L))
        samps <- samps[which(pip==1),]

    }
    samps
}

.check_bisquare_args <- function(manifold,loc,R) {
    stopifnot(is.matrix(loc))
    stopifnot(dimensions(manifold) == ncol(loc))
    stopifnot(is.numeric(R))
    stopifnot(R > 0)
}

# Gaussian Basis Function
.GRBF_wrapper <- function(manifold,mu,std) {
    .check_bisquare_args(manifold,mu,std)
    function(s) {
        stopifnot(ncol(s) == dimensions(manifold))
        dist_sq <- distance(manifold,s,mu)^2
        exp(-0.5* dist_sq/(std^2) )
    }
}

# Bisquare Basis Function
.bisquare_wrapper <- function(manifold,c,R) {
    .check_bisquare_args(manifold,c,R)
    function(s) {
        stopifnot(ncol(s) == dimensions(manifold))
        y <- distance(manifold,s,c)
        (1-(y/R)^2)^2 * (y < R)
    }
}

# Exponential Basis Function
.exp_wrapper <- function(manifold,c,tau) {
    .check_bisquare_args(manifold,c,tau)
    function(s) {
        stopifnot(ncol(s) == dimensions(manifold))
        y <- distance(manifold,s,c)
        exp(-y/tau)
    }
}

# Exponential Basis Function
.Matern32_wrapper <- function(manifold,c,kappa) {
    .check_bisquare_args(manifold,c,kappa)
    function(s) {
        stopifnot(ncol(s) == dimensions(manifold))
        y <- distance(manifold,s,c)
        (1 + sqrt(3)*y/kappa)*exp(-sqrt(3)*y/kappa)
    }
}

#' @rdname concat
#' @aliases concat,Basis-method
#' @noRd
setMethod("concat",signature = "Basis",function(...) {
    l <- list(...)
    if(length(l) < 2)
        stop("Need more than one basis set to concatenate")
    if(!(length(unique(sapply(sapply(l,manifold),type))) == 1))
        stop("Basis need to be on the same manifold")
    G <- l[[1]]

    for (i in 2:length(l)) {
        G@fn <- c(G@fn, l[[i]]@fn)
        G@pars <- c(G@pars, l[[i]]@pars)
        G@df <- rbind(G@df, l[[i]]@df)
    }
    G@n <- length(G@fn)
    G
})

#' @rdname nbasis
#' @aliases nbasis,Basis_obj-method
setMethod("nbasis",signature(.Object="Basis_obj"),function(.Object) {return(.Object@n)})


#' @rdname nbasis
#' @aliases nbasis,SRE-method
setMethod("nbasis",signature(.Object="SRE"),function(.Object) {return(nbasis(.Object@basis))})

#' @aliases count_res,SRE-method
setMethod("count_res",signature="SRE",function(.Object) {
    count_res(.Object@basis)
})

#' @aliases count_res,TensorP_Basis-method
setMethod("count_res",signature="TensorP_Basis",function(.Object) {
    res <- NULL # suppress bindings (it's in the data frame)
    c1 <-  count(.Object@Basis1@df,res)
    c2 <-  count(.Object@Basis2@df,res)

    c_all <- NULL
    max_res_c1 <- c1$res[1] - 1
    for( i in 1:nrow(c2)) {
        new_res <- (max_res_c1 + 1):(max_res_c1 + nrow(c1))
        temp_c1 <- c1
        temp_c1$res <- new_res
        temp_c1$n <- temp_c1$n * c2$n[i]
        c_all <- rbind(c_all,temp_c1)
        max_res_c1 <- max(c_all$res)
    }
    c_all
})

#' @aliases count_res,Basis-method
setMethod("count_res",signature="Basis",function(.Object) {
    res <- NULL # suppress bindings (it's in the data frame)
    count(.Object@df,res)
})


