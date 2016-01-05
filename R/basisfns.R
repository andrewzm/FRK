#' @title Construct a set of radial basis functions
#' @description Construct a set of radial basis functions based on pre-specified location and scale parameters.
#' @param manifold object of class \code{manifold}, for example, a sphere
#' @param loc a matrix of size \code{n} by \code{dimensions(manifold)} indicating centres of basis functions
#' @param scale vector of length \code{n} containing the scale parameters of the basis functions. See details
#' @param type either ``Gaussian'', ``bisquare,'' ``exp,'' or ``Matern32''.
#' @details This functions lays out radial basis functions in a domain of interest based on pre-specified location and scale parameters. If the type is ``Gaussian'', then the scale corresponds to a distance of one standard deviation. If the type is ``bisquare'', then the scale corresponds to the range of support of the bisquare function.
#' @examples
#' library(ggplot2)
#' G <-  radial_basis(manifold = real_line(),
#'                    loc=matrix(1:10,10,1),
#'                    scale=rep(2,10),
#'                    type="bisquare")
#' show_basis(G)
#' @export
radial_basis <- function(manifold=sphere(),loc=matrix(c(1,0),nrow=1),scale=1,type="Gaussian") {
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
#' @description Generate automatically a set of radial basis functions in the domain, and automatically prune in regions of sparse data.
#' @param m object of class \code{manifold}, for example, \code{sphere} or \code{plane}
#' @param data oject of class \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame} containing the data on which basis-function placement is based, see details
#' @param regular an integer indicating the number of regularly-placed basis functions at the first resolution. In two dimensions, this dictates smallest number of basis functions in a row or column at the lowest resolution. If \code{regular=0}, an irregular grid is used, one that is based on the triangulation of the domain with increased mesh density in areas of high data density, see details
#' @param nres the number of basis-function resolutions to use
#' @param prune a threshold parameter which dictates when a basis function is considered irrelevent or unidentifiable, and thus removed, see details
#' @param subsamp the maximum amount of data points to consider when carrying out basis-function placement: these data objects are randomly sampled from the full dataset. Keep this number fairly high (on the order of 10^5) otherwise high resolution basis functions may be spuriously removed
#' @param type the type of basis functions to use. Currently only ``Gaussian'' and ``bisquare'' are supported
#' @details This function automatically places basis functions within the domain of interest. If the domain is a plane or the real line, then the object \code{data} is used to establish the domain boundary.
#'
#' If the manifold is the real line, the basis functions are placed regularly inside the domain, and the number of basis functions at the lowest resolution is dictated by the integer parameter \code{regular} which has to be greater than zero. Each subsequent resolution has twice as many basis functions. The scale of the basis function is set based on the minimum distance between the centre locations following placement. The scale is equal to the minimum distance if the type of basis function is Gaussian, and 1.5x this value if the function is bisquare. See \code{\link{radial_basis}} for details.
#'
#' If the manifold is a plane, and \code{regular > 0}, then basis functions are placed regularly within the bounding box of \code{data}, with the smallest number of basis functions in each row or column equal to the value of \code{regular} in the lowest resolution. Subsequent resolutions have twice the number of basis functions in each row or column. If \code{regular = 0}, then the function \code{INLA::inla.nonconvex.hull} is used to construct a (non-convex) hull around the data. The buffer and smoothness of the hull is determined by the parameter \code{convex}. Once the domain boundary is found,  \code{INLA::inla.mesh.2d} is used to construct a triangular mesh such that the node vertices coincide with data locations, subject to some minimum and maximum triangular side length constraints. The result is a mesh which is dense in regions of high data density and not dense in regions of sparse data. Even in this case, the scale is taken to be the minimum distance between basis function centres. This may be changed in a future revision.
#'
#' If the manifold is a sphere, then basis functions are placed on the centroids of the discrete global grid (DGG), with the first basis resolution corresponding to the second resolution of the DGG (32 globally).  It is not recommended to go above \code{nres == 4} which contains 812 locations (for a total of 1208 basis functions).
#'
#' Basis functions that are not influenced by data points may hinder convergence of the EM algorithm, since the associated hidden states are by and large unidentifiable. We hence provide a means to automatically remove such basis functions through the parameter \code{prune}. The final set only contains basis functions for which the column sums in the associated matrix \eqn{S} (which, recall, is the value/average of the basis functions at/over the data points/polygons) is greater than \code{prune}. If \code{prune == 0}, no basis functions are removed from the original design.
#' @examples
#' ### Load and process the AIRS dataset
#'
#' library(dplyr)
#' library(sp)
#' library(ggplot2)
#'
#' ### Create a synthetic dataset
#' d <- data.frame(lon = runif(n=1000,min = -179, max = 179),
#'                 lat = runif(n=1000,min = -90, max = 90),
#'                 z = rnorm(5000))
#' coordinates(d) <- ~lon + lat
#' proj4string(d)=CRS("+proj=longlat")
#'
#' ### Now create basis functions over sphere
#' G <- auto_basis(m = sphere(),data=d,
#'                 nres = 2,prune=15,
#'                 type = "bisquare",
#'                 subsamp = 20000)
#'
#' ### Plot and note how some basis functions are removed in Antarctica
#' show_basis(G,draw_world())
#'
#' @export
auto_basis <- function(m = plane(),data,regular=1,nres=2,prune=0,subsamp=10000,type="Gaussian") {
    if(!is(m,"manifold")) stop("m needs to be a manifold")
    if(!is.numeric(nres) | nres < 0) stop("nres needs to be greater than zero")
    if(!is.numeric(prune) | prune < 0) stop("prune needs to be greater than zero")
    if(!is.numeric(subsamp) | subsamp < 0) stop("subsamp needs to be greater than zero")
    if(!type %in% c("bisquare","Gaussian","exp","Matern32")) stop("type of basis functions must be 'Gaussian,' 'bisquare,' 'exp,' or 'Matern32'.")
    if((is(m,"sphere")  | is(m,"real_line")) & regular == 0) stop("Irregular basis only available on planes")

    isea3h <- centroid <- res <- NULL #(suppress warnings, these are loaded from data)
    coords <- coordinates(data)

    if(is(m,"plane") & regular==0) {
        if(!requireNamespace("INLA")) stop("For automatic basis generation INLA needs to be installed for constructing basis function centres. Please install it using install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/stable\")")
    }

    if(nrow(coords)>subsamp) {
        coords <- coords[sample(nrow(coords),size=subsamp,replace=FALSE),]
    }

    xrange <- range(coords[,1])
    yrange <- range(coords[,2])

    if(is(m,"plane") & regular == 0) bndary_seg = INLA::inla.nonconvex.hull(coords,convex=-0.05)
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

    #if(is(m,"sphere")) load(system.file("extdata","isea3h.rda", package = "FRK"))
    if(is(m,"sphere")) {
        isea3h <- load_dggrids(res = nres)
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
            ## Generate mesh and use these as centres
            xgrid <- seq(xrange[1] + diff(xrange)/(2*nx*i), xrange[2] - diff(xrange)/(2*nx*i), length =nx*i)
            ygrid <- seq(yrange[1] + diff(yrange)/(2*ny*i), yrange[2] - diff(yrange)/(2*ny*i), length =ny*i)

            this_res_locs <- xgrid %>%
                expand.grid(ygrid) %>%
                as.matrix()
        } else if(is(m,"real_line")) {
            this_res_locs <- matrix(seq(xrange[1],xrange[2],length=i*regular))
        } else if(is(m,"sphere")) {
            this_res_locs <- as.matrix(filter(isea3h,centroid==1,res==(i-1))[c("lon","lat")])
        }
        ## Set scales: To 1.5x the distance to nearest basis
        ## Refine: Remove basis which are not influenced by data and re-find the scales
        for(j in 1:2) {
            D <- FRK::distance(m,this_res_locs,this_res_locs)
            if(nrow(D) == 1) {
                this_res_scales <-max(diff(xrange),diff(yrange))/2
            } else {
                diag(D) <- Inf
                this_res_scales <- apply(D,1,min)
            }

            this_res_basis <- radial_basis(manifold = m,
                                           loc=this_res_locs,
                                           scale=ifelse(type=="Gaussian",1,1.5)*this_res_scales,
                                           type=type)
            if(prune > 0) {
                if(j==1) {
                    rm_idx <- which(colSums(eval_basis(this_res_basis,coords)) < prune)
                    if(length(rm_idx) == length(this_res_scales)) stop("prune is too large -- all functions at a resolution removed. Consider also removing number of resolutions.")
                    if(length(rm_idx) >0) this_res_locs <- this_res_locs[-rm_idx,,drop=FALSE]
                }
            } else {
                break
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
    if(G_basis@n > 2000) warning("More than 2000 basis functions")
    G_basis
}

#' @title Add the time coordinate to 2D spatial basis functions
#' @description Given a set of 2D spatial basis functions and a vector of knots in time, this function repeats the spatial basis at every temporal knot, adding the third dimension to the centroid as appopriate.
#' @param G_spatial an object of class Basis on a 2D manifold
#' @param t_knots a vector of numbers locating the knots in time
#' @param manifold a 3D space-time manifold, typically STsphere or STplane
#' @examples
#' G_spatial <-  radial_basis(manifold = sphere(),
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
    new("TensorP_Basis",Basis1=Basis1, Basis2=Basis2, n = Basis1@n * Basis2@n)
})


#' @rdname eval_basis
#' @aliases eval_basis,Basis-matrix-method
setMethod("eval_basis",signature(basis="Basis",s="matrix"),function(basis,s,output = "matrix"){
    stopifnot(output %in% c("list","matrix"))
    space_dim <- dimensions(manifold(basis))
    .point_eval_fn(basis@fn,s[,1:space_dim,drop=F],output)
})

#' @rdname eval_basis
#' @aliases eval_basis,Basis-SpatialPointsDataFrame-method
setMethod("eval_basis",signature(basis="Basis",s="SpatialPointsDataFrame"),function(basis,s,output = "matrix"){
    stopifnot(output %in% c("matrix","list"))
    space_dim <- dimensions(manifold(basis))
    .point_eval_fn(basis@fn,coordinates(s)[,1:space_dim,drop=F],output)
})

#' @rdname eval_basis
#' @aliases eval_basis,Basis-SpatialPolygonsDataFrame-method
setMethod("eval_basis",signature(basis="Basis",s="SpatialPolygonsDataFrame"),function(basis,s,output = "matrix"){
    stopifnot(output %in% c("matrix","list"))
    X <- list()
    print("Averaging over polygons")

    if(opts_FRK$get("parallel") > 0) {
        cl <- makeCluster(opts_FRK$parallel)
        X <- mclapply(1:length(s), function(i) {
            samps <- .samps_in_polygon(basis,s,i)
            colSums(.point_eval_fn(basis@fn,samps))/nrow(samps)
        },mc.cores = opts_FRK$get("parallel"))
        stopCluster(cl)
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
    i <- 1 #suppress binding warning

    S <- list()
    for(i in 1:ncol(S1)) {
        S[[i]] <- S1[,i] * S2
    }
    S <- do.call("cBind",S)
    S <- as(S,"dgCMatrix")

    # S <- foreach(i = 1:ncol(S1),.combine="cBind") %do% {
    #     S1[,i] * S2
    # } %>% as("dgCMatrix")

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

    i <- 1 #suppress binding warning

    S <- list()
    for(i in 1:ncol(S1)) {
        S[[i]] <- S1[,i] * S2
    }
    S <- do.call("cBind",S)
    S <- as(S,"dgCMatrix")

    # S <- foreach(i = 1:ncol(S1),.combine="cBind") %do% {
    #     S1[,i] * S2
    # } %>% as("dgCMatrix")

    S
})



#' @rdname eval_basis
#' @aliases eval_basis,Basis-STIDF-method
setMethod("eval_basis",signature(basis="Basis",s="STIDF"),function(basis,s,output = "matrix"){
    stopifnot(output %in% c("matrix","list"))
    space_dim <- dimensions(manifold(basis))
    .point_eval_fn(basis@fn,cbind(coordinates(s),t=s@data$t)[,1:space_dim,drop=F],output)
})






.point_eval_fn <- function(flist,s,output="matrix") {
    #if(!opts_FRK$get("Rhipe")) {
    if(1) {
        #x <- sapply(flist,function(f) f(s))
        x <- do.call("cbind",sapply(flist,function(f) f(s),simplify=FALSE))
        as(x,"Matrix")
    } else { ## The below works but likely to be slower.. whole prediction should be parallelised
        envlist <- lapply(flist, function(f) environment(f))
        flist <- lapply(flist, function(f) parse(text = deparse(f)))
        x <- rhwrapper(Ntot = nrow(s),
                       N = 4000,
                       type="Matrix",
                       f_expr = .rhpoint_eval_fn,
                       flist=flist,
                       envlist = envlist,
                       s=s)
        as(data.matrix(x),"Matrix")
    }
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
        #             ggplot(X) + geom_point(aes(X1,X2),pch=".") + geom_path(data=coords,aes(lon,lat)) + coord_map("ortho",xlim=c(-180,180),ylim=c(-90,90),orientation=c(-180,-45,45))
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


