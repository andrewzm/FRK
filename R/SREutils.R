#' @title Construct SRE object
#' @export
SRE <- function(f,data,basis,BAUs,est_error=T) {

    .check_args(f=f,data=data,basis=basis,BAUs=BAUs,est_error=est_error)
    av_var <-all.vars(f)[1]
    ndata <- length(data)

    S <- Ve <- Vfs <- X <- Z <- Cmat <- list()
    for(i in 1:ndata) {
        data_proc <- map_data_to_BAUs(data[[i]],
                                 BAUs,
                                 av_var = av_var,
                                 variogram.formula = f,
                                 est_error=est_error)

        L <- .gstat.formula(f,data=data_proc)
        X[[i]] <- as(L$X,"Matrix")
        Z[[i]] <- Matrix(L$y)
        Ve[[i]] <- Diagonal(x=data_proc$std^2)

        if(is(data_proc,"SpatialPolygonsDataFrame")) {
            data_proc$id <- 1:length(data_proc)
            overlap <- over(SpatialPoints(coordinates(BAUs)),data_proc)
            i_idx <- as.numeric(na.exclude(overlap$id))
            j_idx <- which(!is.na(overlap$id))
        } else {
             i_idx <- 1:length(data_proc)
             j_idx <- which(row.names(BAUs) %in% row.names(data_proc))
        }
        Cmat[[i]] <- sparseMatrix(i=i_idx,j=j_idx,x=1,dims=c(length(data_proc),nrow(BAUs)))
        Cmat[[i]] <- Cmat[[i]] / rowSums(Cmat[[i]])
        Vfs[[i]] <- Diagonal(x=as.numeric(Cmat[[i]]^2 %*% BAUs$fs )) # Assuming no obeservations overlap
        S[[i]] <- eval_basis(basis, s = data_proc)

    }

    S <- do.call("rBind",S)
    X <- do.call("rBind",X)
    Z <- do.call("rBind",Z)
    Ve <- do.call("bdiag",Ve)
    Vfs <- do.call("bdiag",Vfs)

    new("SRE",
        data=data,
        basis=basis,
        BAUs=BAUs,
        f = f,
        S = S,
        Ve = Ve,
        Vfs = Vfs,
        X = X,
        Z = Z,
        Cmat = do.call("rBind",Cmat),
        mu_eta = Matrix(0,nbasis(basis),1),
        S_eta = Diagonal(x = rep(1,nbasis(basis))),
        alphahat = solve(t(X) %*% X) %*% t(X) %*% Z,
        Khat = Diagonal(n=nbasis(basis),x = var(Z[,1])),
        sigma2fshat = mean(diag(Ve)) / mean(diag(Vfs)))
}

#' @title Estimate SRE model parameters
#' @export
SRE.fit <- function(SRE_model,n_EM = 100L, tol = 1e-5, method="EM",print_lik=F) {
    if(!is.numeric(n_EM)) stop("n_EM needs to be an integer")
    if(!(n_EM <- round(n_EM)) > 0) stop("n_EM needs to be greater than 0")
    if(!is.numeric(tol)) stop("tol needs to be a number greater than zero")
    if(!(tol > 0)) stop("tol needs to be a number greater than zero")
    if(!(method == "EM")) stop("Currently only the EM algorithm is implemented for parameter estimation")
    if(!(is.logical(print_lik))) stop("print_lik needs to be a logical quantity")

    n <- nbasis(SRE_model)
    X <- SRE_model@X
    lk <- rep(0,n_EM)


    if(opts_FRK$get("progress")) pb <- txtProgressBar(min = 0, max = n_EM, style = 3)
    for(i in 1:n_EM) {
        lk[i] <- .loglik(SRE_model)  # Compute likelihood as in Katzfuss and Cressie (2011)
        SRE_model <- .SRE.Estep(SRE_model)
        SRE_model <- .SRE.Mstep(SRE_model)
        if(opts_FRK$get("progress")) setTxtProgressBar(pb, i)
        if(i>2) if(lk[i] - lk[i-1] < tol) {
            print("Minimum tolerance reached")
            break
        }
    }
    if(opts_FRK$get("progress")) close(pb)

    if(i == n_EM) print("Maximum EM iterations reached")
    if(print_lik) {
        plot(lk[1:i][-1],ylab="log likelihood")
        print("Warning: Ignoring constants in log-likelihood computation")
    }
    SRE_model
}

#setMethod("SRE.predict",signature(Sm="SRE", pred_locs="NULL",depname="character"),

#' @title Predict using SRE model
#' @export
SRE.predict <- function(Sm,pred_locs = Sm@BAUs,use_centroid=T) {
### CHANGE INPUT TO IDX!!! OR NAMES OR SOMETHING SO THAT WE CAN PREDICT ON A SUBSET OF BAUs!!
    if(!(is(pred_locs,"SpatialPolygonsDataFrame"))) stop("Predictions need to be over BAUs or spatial polygons")
    if(!("fs" %in% names(pred_locs))) {
        warning("data should contain a field 'fs' containing a basis function for fine-scale variation. Setting basis function equal to one everywhere.")
        pred_locs$fs <- 1
    }
    if(!(all(pred_locs$fs > 0))) stop("fine-scale variation basis function needs to be nonnegative everywhere")

    if(opts_FRK$get("Rhipe")) {
        pred_locs <- rhwrapper(Ntot = length(pred_locs),
                               N = 4000,
                               f_expr = .rhSRE.predict,
                               Sm = Sm,
                               pred_locs = pred_locs,
                               use_centroid = use_centroid)
    } else {

        pred_locs <- .SRE.predict(Sm=Sm,pred_locs=pred_locs,use_centroid=use_centroid)
    }
    pred_locs
}


.SRE.predict <- function(Sm,pred_locs,use_centroid) {
        depname <- all.vars(Sm@f)[1]
        pred_locs[[depname]] <- 0.1
        L <- .gstat.formula(Sm@f,data=pred_locs)
        X = as(L$X,"Matrix")
        if(use_centroid) {
            S0 <- eval_basis(Sm@basis,as.matrix(pred_locs[coordnames(Sm@data[[1]])]@data))
        } else {
            S0 <- eval_basis(Sm@basis,pred_locs)
        }
        pred_locs[[depname]] <- NULL

        alpha <- Sm@alphahat
        K <- Sm@Khat
        sigma2fs <- Sm@sigma2fshat
        D <- sigma2fs*Sm@Vfs + Sm@Ve
        Dinv <- Diagonal(x=1/diag(D))

        idx <- match(row.names(pred_locs),row.names(Sm@BAUs))
        sig2_Vfs_pred <- Diagonal(x=sigma2fs*pred_locs$fs)

        idx <- match(row.names(pred_locs),row.names(Sm@BAUs))
        C <- Sm@Cmat[,idx]

        LAMBDA <- as(bdiag(Sm@Khat,sig2_Vfs_pred),"symmetricMatrix")
        LAMBDAinv <- chol2inv(chol(LAMBDA))
        PI <- cBind(S0, .symDiagonal(n=nrow(pred_locs)))
        Qx <- t(PI) %*% t(C) %*% solve(Sm@Ve) %*% C %*% PI + LAMBDAinv
        temp <- cholPermute(Qx)
        ybar <- t(PI) %*%t(C) %*% solve(Sm@Ve) %*% (Sm@Z - C %*% X %*% alpha)
        x_mean <- cholsolve(Qx,ybar,perm=T,cholQp = temp$Qpermchol, P = temp$P)
        Partial_Cov <- Takahashi_Davis(Qx,cholQp = temp$Qpermchol,P = temp$P)
        x_margvar <- diag(Partial_Cov)

        pred_locs[["mu"]] <- as.numeric(X %*% alpha + PI %*% x_mean)

        ## variance to hard to compute all at once -- do it in blocks of 1000
        temp <- rep(0,length(pred_locs))
        batching=cut(1:nrow(PI),breaks = seq(0,nrow(PI)+1000,by=1000),labels=F)
        for(i in 1:max(unique(batching))) {
            idx = which(batching==i)
            temp[idx] <- as.numeric(rowSums((PI[idx,] %*% Partial_Cov)*PI[idx,]))
        }
        pred_locs[["var"]] <- temp
        #pred_locs[["var"]] <- as.numeric(rowSums((PI %*% Partial_Cov)*PI))
        pred_locs
}

.SRE.Estep <- function(Sm) {

    alpha <- Sm@alphahat
    K <- Sm@Khat
    sigma2fs <- Sm@sigma2fshat

    D <- sigma2fs*Sm@Vfs + Sm@Ve
    Dinv <- Diagonal(x=1/diag(D))

    Kinv <- chol2inv(chol(K))

    S_eta <- chol2inv(chol(crossprod(sqrt(Dinv) %*% Sm@S) + Kinv))
    mu_eta <- S_eta %*% (t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha))

    Sm@mu_eta <- mu_eta
    Sm@S_eta <- S_eta

    Sm
}


.SRE.Mstep <- function(Sm) {

    mu_eta <- Sm@mu_eta
    S_eta <- Sm@S_eta
    alpha_init <- Sm@alphahat
    sigma2fs_init <- Sm@sigma2fshat

    K <- S_eta + tcrossprod(mu_eta)
    alpha <- alpha_init
    sigma2fs <- sigma2fs_init
    converged <- FALSE

    Omega_diag1 <- diag2(Sm@S %*% as.matrix(S_eta),t(Sm@S)) +
                   diag2(Sm@S %*% mu_eta %*% t(mu_eta), t(Sm@S))

    while(!converged) {
        J <- function(sigma2fs) {
            if(sigma2fs < 0) {
                return(Inf)
            } else {
                D <- sigma2fs*Sm@Vfs + Sm@Ve
                Dinv <- solve(D)

                -(-0.5*tr(Dinv %*% Sm@Vfs) +
                      0.5*tr(Dinv %*% Sm@Vfs %*% Dinv %*% Omega_diag)
                )
            }
        }

        resid <- Sm@Z - Sm@X %*% alpha
        Omega_diag <- Omega_diag1 -
            2*diag2(Sm@S %*% mu_eta, t(resid)) +
            diag2(resid,t(resid))
        Omega_diag <- Diagonal(x=Omega_diag)

        # Repeat until finding values on opposite sides of zero
        amp_factor <- 10; OK <- 0
        while(!OK) {
            amp_factor <- amp_factor * 10
            if(!(sign(J(sigma2fs/amp_factor)) == sign(J(sigma2fs*amp_factor)))) OK <- 1
        }
        sigma2fs_new <- uniroot(f = J,interval = c(sigma2fs/amp_factor,sigma2fs*amp_factor))$root
        D <- sigma2fs_new*Sm@Vfs + Sm@Ve
        Dinv <- Diagonal(x=1/diag(D))
        alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*% t(Sm@X) %*% Dinv %*% (Sm@Z - Sm@S %*% mu_eta)

        if(max(sigma2fs_new / sigma2fs, sigma2fs / sigma2fs_new) < 1.001) converged <- TRUE
        sigma2fs <- sigma2fs_new
    }

    Sm@Khat <- K
    Sm@alphahat <- alpha
    Sm@sigma2fshat <- sigma2fs

    Sm
}





.loglik <- function(Sm) {

    D <- Sm@sigma2fshat*Sm@Vfs + Sm@Ve
    Dinv <- solve(D)
    r <- Sm@Z - Sm@X %*% Sm@alphahat  - Sm@S %*% Sm@mu_eta

    as.numeric(-0.5*determinant(D)$modulus -
                   0.5 * t(r) %*% Dinv %*% r)
}


.check_args <- function(f,data,basis,BAUs,est_error) {
    if(!is(f,"formula")) stop("f needs to be a formula.")
    if(!all(all.vars(f)[-1] %in% names(BAUs))) stop("All covariates need to be in the SpatialPolygons BAU object")
    if(!is(data,"list")) stop("Please supply a list of Spatial objects.")
    if(!all(sapply(data,function(x) is(x,"Spatial")))) stop("All data list elements need to be of class Spatial")
    if(!all(sapply(data,function(x) all.vars(f)[1] %in% names(x)))) stop("All data list elements to have values for the dependent variable")
    if(!all(sapply(data,function(x) identical(proj4string(x), proj4string(BAUs))))) stop("Please ensure all data items and BAUs have the same coordinate reference system")
    if(!is(basis,"Basis")) stop("basis needs to be of class Basis (package FRK)")
    #if(is(data,"SpatialPolygonsDataFrame") & !("std" %in% names(data))) stop("Polygon data needs to contain a field 'std' denoting the observation error")
    if(!("fs" %in% names(BAUs))) {
        warning("BAUs should contain a field 'fs' containing a basis function for fine-scale variation. Setting basis function equal to one everywhere.")
        BAUs$fs <- 1
    }
    if(!is(BAUs,"SpatialPolygonsDataFrame")) stop("fine-scale variation basis function needs to be nonnegative everywhere")
    if((is(basis@manifold,"sphere")) & !all((coordnames(BAUs) == c("lon","lat")))) stop("Since a sphere is being used, please ensure that all coordinates (including those of BAUs) are in (lon,lat)")
    if(!est_error & !all(sapply(data,function(x) "std" %in% names(x)))) stop("If observational error is not going to be estimated, please supply a field 'std' in the data objects")
}

.gstat.formula <- function (formula, data)
{
    if (is(data, "SpatialPixels") && anyDuplicated(data@grid.index) !=
        0)
        gridded(data) = FALSE
    m = model.frame(terms(formula), as(data, "data.frame"), na.action = na.fail)
    Y = model.extract(m, "response")
    if (length(Y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)
    has.intercept = attr(Terms, "intercept")
    if (gridded(data))
        grid = gridparameters(data)
    else grid = numeric(0)
    xlevels = .getXlevels(Terms, m)
    list(y = Y, locations = coordinates(data), X = X, call = call,
         has.intercept = has.intercept, grid = as.double(unlist(grid)),
         xlevels = xlevels)
}
