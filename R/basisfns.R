#' @title Construct a set of radial basis functions
#' @export
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
            fn[[i]] <-  .GRBF_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        } else {
            fn[[i]] <-  .bisquare_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        }
        pars[[i]] <- list(loc = matrix(loc[i,],nrow=1), scale=scale[i])
    }
    df <- data.frame(loc,scale,res=1)
    this_basis <- new("Basis", manifold=manifold, pars=pars, n=n, fn=fn, df=df)
    return(this_basis)
}

#' @title Automatic basis-function placement
#' @export
auto_basis <- function(m = plane(),data,regular=1,nres=2,prune=0,subsamp=10000,type="Gaussian") {
    if(!is(m,"manifold")) stop("m needs to be a manifold")
    if(!is.numeric(nres) | nres < 0) stop("nres needs to be greater than zero")
    if(!is.numeric(prune) | prune < 0) stop("prune needs to be greater than zero")
    if(!is.numeric(subsamp) | subsamp < 0) stop("subsamp needs to be greater than zero")
    if(!type %in% c("bisquare","Gaussian")) stop("type of basis functions can only be Gaussian or bisquare")
    if((is(m,"sphere")  | is(m,"real_line")) & regular == 0) stop("Irregular basis only available on planes")

    isea3h <- centroid <- res <- NULL #(suppress warnings, these are loaded from data)
    coords <- coordinates(data)

    if(is(m,"plane")) {
        if(!require(INLA)) stop("For automatic basis generation INLA needs to be installed. Please install it using install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/stable\")")
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
    if(is(m,"sphere")) data(isea3h, envir=environment())

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
    G_basis
}


#' @rdname eval_basis
#' @aliases eval_basis,Basis-matrix-method
setMethod("eval_basis",signature(basis="Basis",s="matrix"),function(basis,s,output = "matrix"){
    stopifnot(output %in% c("list","matrix"))
    space_dim <- dimensions(manifold(basis))
    .point_eval_fn(basis@fn,s[,1:space_dim],output)
})

#' @rdname eval_basis
#' @aliases eval_basis,Basis-SpatialPointsDataFrame-method
setMethod("eval_basis",signature(basis="Basis",s="SpatialPointsDataFrame"),function(basis,s,output = "matrix"){
    stopifnot(output %in% c("matrix","list"))
    space_dim <- dimensions(manifold(basis))
    .point_eval_fn(basis@fn,coordinates(s)[,1:space_dim],output)
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


.point_eval_fn <- function(flist,s,output="matrix") {
    #if(!opts_FRK$get("Rhipe")) {
    if(1) {
        x <- sapply(flist,function(f) f(s))
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
        #             ggplot(X) + geom_point(aes(X1,X2),pch=".") + coord_map() + geom_path(data=coords,aes(lon,lat)) + coord_map("ortho",xlim=c(-180,180),ylim=c(-90,90),orientation=c(-180,-45,45))
    }
    samps
}


.GRBF_wrapper <- function(manifold,mu,std) {
    stopifnot(is.matrix(mu))
    stopifnot(dimensions(manifold) == ncol(mu))
    stopifnot(is.numeric(std))
    stopifnot(std > 0)
    function(s) {
        stopifnot(ncol(s) == dimensions(manifold))
        dist_sq <- distance(manifold,s,mu)^2
        exp(-0.5* dist_sq/(std^2) )
    }
}

# Gaussian Asymmetric Basis Function
.bisquare_wrapper <- function(manifold,c,R) {
    stopifnot(dimensions(manifold) == ncol(c))
    stopifnot(is.numeric(R))
    stopifnot(R > 0)
    function(s) {
        y <- distance(manifold,s,c)
        (1-(y/R)^2)^2 * (y < R)
    }
}



#' @rdname concat
#' @aliases concat,Basis-method
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
#' @aliases nbasis,Basis-method
setMethod("nbasis",signature(.Object="Basis"),function(.Object) {return(.Object@n)})

#' @rdname nbasis
#' @aliases nbasis,SRE-method
setMethod("nbasis",signature(.Object="SRE"),function(.Object) {return(nbasis(.Object@basis))})


