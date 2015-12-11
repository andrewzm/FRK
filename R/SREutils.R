#' @title Construct SRE object, fit and predict
#' @description Main constructor of spatial random effects (SRE) object. Please see \code{\link{SRE-class}} for more details on the object's properties and methods.
#' @param f \code{R} formula relating the dependent variable (or transformations thereof) to covariates
#' @param data list of objects of class \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}
#' @param basis object of class \code{Basis}
#' @param BAUs object of class \code{SpatialPolygonsDataFrame}, the data frame which must contain covariate information as well as a field \code{fs} describing the fine-scale variation up to a constant of proportionality
#' @param est_error flag indicating whether the variance of the data should be estimated from variogram techniques. If this is set to 0, then \code{data} must contain a field \code{std}
#' @param SRE_model object returned from the constructor \code{SRE()}
#' @param n_EM maximum number of iterations for the EM algorithm
#' @param tol convergence tolerance for the EM algorithm
#' @param method parameter estimation method to employ. Currently only \code{EM} is supported
#' @param print_lik flag indicating whether likelihood should be printed or not on convergence of the estimation algorithm
#' @param pred_locs object of class \code{SpatialPolygonsDataFrame} containing the prediction polygons. The data frame of \code{pred_locs} must contain covariate information as well as a field \code{fs} describing the fine-scale variation up to a constant of proportionality. Ideally, \code{pred_locs} is identical to \code{BAUs} above
#' @param use_centroid a flag indicating whether the prediction over a BAU can be simply taken as a point prediction at the BAU's centroid. This should only be done if the BAUs on which the model is trained coincide with the BAUs used in \code{SRE()}
#' @details \code{SRE()} is the main function in the program as it constructs a spatial random effects model from the user-defined formula, data object, basis functions and a set of basic aerial units (BAUs). The function first takes each object in the list \code{data} and maps it to the BAUs -- this entails binning the point-referenced data into BAUs (and averaging within the BAU) and finding which BAUs are influenced by the polygon datasets. Following this the incidence matrix \code{Cmat} is constructed, which appears in the observation model \eqn{Z = CY + e}, where \eqn{C} is the incidence matrix.
#'
#' The SRE model is given by \eqn{Y = X\beta + S\eta + v} where \eqn{X} are the covariates at BAU level, $\beta$ are the regression coefficients, \eqn{S} are the basis functions evaluated at the BAU level, \eqn{\eta} are the basis function weights, and \eqn{v} is the fine scale variation (at the BAU level). The covariance matrix of \eqn{v} is diagonal and proportional to the field 'fs' in the BAUs (typically set to one). The constant of proportionality is estimated in the EM algorithm. All required matrices (\eqn{S,X} etc.) are computed and returned as part of the object, please see \code{\link{SRE-class}} for more details.
#'
#'\code{SRE.fit()} takes an object of class \code{SRE} and estimates all unknown parameters, namely the covariance matrix \eqn{K}, the fine scale variance \eqn{\sigma^2_{fs}} and the regression parameters \eqn{\alpha}. The only method currently implemented is the expectation maximisation (EM) algorithm, which the user configures through \code{n_EM} and \code{tol}. The latter parameter, \code{tol}, is used as in Katzfuss and Cressie (2011), that is, the log-likelihood (given in Equation (16) in that work) is evaluated at each iteration at the current parameter estimate and convergence is assumed reach when this quantity stops changing by more than \code{tol}.
#'
#'The actual computations for the E-step and M-step are relatively straightforward. The E-step contains an inverse of an \eqn{n \times n} matrix, where \code{n} is the number of basis functions which should not exceed 2000. The M-step first updates the matrix \eqn{K}, which only depends on the sufficient statistics of the basis weights \eqn{\eta}. Then, a line search is used to update the fine-scale variance \eqn{\sigma^2_{fs}} and the regression parameters \eqn{\alpha}. Since the udpates of these last two parameters are not independent, the updates are iterated until the change in \eqn{\sigma^2_{fs}} is no more than 0.1\%.
#'
#'Once the parameters are fitted, the \code{SRE} object is passed onto the function \code{SRE.predict()} in order to carry out optimal predictions over selected polygons, typically the same BAUs used to construct the SRE model with \code{SRE()}. The first part of the prediction process is to construct the matrix \eqn{S} by averaging the basis functions over the prediction polygons, using Monte Carlo integration with 1000 samples. This is a computationally-intensive process and can be distributed over the Hadoop backend if desired. On the other hand, one could set \code{use_centroid = TRUE} and treat the prediction over the entire polygon as that of a BAU at the centre of the polyon. This will yield valid results only if the polygon is itself a BAU (which many times it is) and if the BAUs are themselves relatively small. Once the matrix \eqn{S} is found, a standard Gaussian inversion using the estimated parameters is used. In order to take advantage of the sparsity of the matrices we carry out sparse matrix inversions using the theory of Erisman and Tinney (1975). For a proof of the validity of this technique please see the accompanying paper.
#'
#'\code{SRE.predict} returns the BAUs, which are of class \code{SpatialPolygonsDataFrame}, with two added attributes, \code{mu} and \code{var}. These can then be easily plotted using \code{spplot} or \code{ggplot2} (in conjunction with \code{\link{SpatialPolygonsDataFrame_to_df}}) as shown in the package vignettes.
#' @references
#' Katzfuss, M., & Cressie, N. (2011). Spatio-temporal smoothing and EM estimation for massive remote-sensing data sets. Journal of Time Series Analysis, 32(4), 430--446.
#' Erisman, A. M., & Tinney, W. F. (1975). On computing certain elements of the inverse of a sparse matrix. Communications of the ACM, 18(3), 177--179.
#' @export
#' @examples
#' library(sp)
#' library(ggplot2)
#' library(dplyr)
#'
#' ### Generate process and data
#' sim_process <- data.frame(x = seq(0.005,0.995,by=0.01)) %>%
#'     mutate(y=0,proc = sin(x*10) + 0.3*rnorm(length(x)))
#' sim_data <- sample_n(sim_process,50) %>%
#'     mutate(z = proc + 0.1*rnorm(length(x)), std = 0.1)
#' coordinates(sim_data) = ~x + y# change into an sp object
#' grid_BAUs <- auto_BAUs(manifold=real_line(),data=sim_data,cellsize = c(0.01),type="grid")
#' grid_BAUs$fs = 1
#'
#' ### Set up SRE model
#' G <- auto_basis(m = real_line(),
#'                 data=sim_data,
#'                 nres = 2,
#'                 regular = 6,
#'                 type = "bisquare",
#'                 subsamp = 20000)
#' f <- z ~ 1
#' S <- SRE(f,list(sim_data),G,
#'          grid_BAUs,
#'          est_error = FALSE)
#'
#' ### Fit with 5 EM iterations so as not to take too much time
#' S <- SRE.fit(S,n_EM = 5,tol = 1e-5,print_lik=TRUE)
#'
#' ### Predict over BAUs
#' grid_BAUs <- SRE.predict(S,pred_locs = grid_BAUs,use_centroid = TRUE)
#'
#' ### Plot
#' X <- slot(grid_BAUs,"data") %>%
#'      filter(x >= 0 & x <= 1)
#' g1 <- LinePlotTheme() +
#'    geom_line(data=X,aes(x,y=mu)) +
#'    geom_errorbar(data=X,aes(x=x,ymax = mu + 2*sqrt(var), ymin= mu - 2*sqrt(var))) +
#'    geom_point(data = data.frame(sim_data),aes(x=x,y=z),size=3) +
#'    geom_line(data=sim_process,aes(x=x,y=proc),col="red")
#' print(g1)
SRE <- function(f,data,basis,BAUs,est_error=TRUE) {

    .check_args(f=f,data=data,basis=basis,BAUs=BAUs,est_error=est_error)
    av_var <-all.vars(f)[1]
    ndata <- length(data)

    S <- Ve <- Vfs <- X <- Z <- Cmat <- list()
    for(i in 1:ndata) {
        if(est_error) data[[i]]$std <- 0 ## Just set it to something, this will be overwritten later on
        data_proc <- map_data_to_BAUs(data[[i]],
                                     BAUs,
                                     av_var = av_var)

        if(any(is.na(data_proc@data)))
            stop("NAs found when mapping data to BAUs. Are you sure all your data are covered by BAUs?")

        if(est_error) {
            if(is(data_proc,"ST")) stop("Estimation of error not yet implemented for spatio-temporal data")
            data_proc <- est_obs_error(data_proc,variogram.formula=f)
        }


        L <- .gstat.formula(f,data=data_proc)
        X[[i]] <- as(L$X,"Matrix")
        Z[[i]] <- Matrix(L$y)
        Ve[[i]] <- Diagonal(x=data_proc$std^2)


        C_idx <- BuildC(data_proc,BAUs)

        Cmat[[i]] <- sparseMatrix(i=C_idx$i_idx,
                                  j=C_idx$j_idx,
                                  x=1,dims=c(length(data_proc),
                                             length(BAUs)))

        Cmat[[i]] <- Cmat[[i]] / rowSums(Cmat[[i]]) ## Average BAUs for polygon observations

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

#' @rdname SRE
#' @export
SRE.fit <- function(SRE_model,n_EM = 100L, tol = 1e-5, method="EM",print_lik=FALSE) {
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
        plot(lk[1:i],ylab="log likelihood")
        print("Warning: Ignoring constants in log-likelihood computation")
    }
    SRE_model
}

#setMethod("SRE.predict",signature(SRE_model="SRE", pred_locs="NULL",depname="character"),

#' @rdname SRE
#' @export
SRE.predict <- function(SRE_model,pred_locs = SRE_model@BAUs,use_centroid=TRUE) {
### CHANGE INPUT TO IDX!!! OR NAMES OR SOMETHING SO THAT WE CAN PREDICT ON A SUBSET OF BAUs!!
    if(is(pred_locs,"Spatial") & !(is(pred_locs,"SpatialPolygonsDataFrame")))
        stop("Predictions need to be over BAUs or spatial polygons")
    if(is(pred_locs,"ST") & !(is(pred_locs,"STFDF")))
       if(!(is(pred_locs@sp,"SpatialPolygonsDataFrame")))
        stop("Predictions need to be over BAUs or STFDFs with spatial polygons")
    if(!("fs" %in% names(pred_locs@data))) {
        warning("data should contain a field 'fs' containing a basis function for fine-scale variation. Setting basis function equal to one everywhere.")
        pred_locs$fs <- 1
    }
    if(!(all(pred_locs$fs > 0))) stop("fine-scale variation basis function needs to be nonnegative everywhere")

    if(opts_FRK$get("Rhipe")) {
        pred_locs <- rhwrapper(Ntot = length(pred_locs),
                               N = 4000,
                               f_expr = .rhSRE.predict,
                               Sm = SRE_model,
                               pred_locs = pred_locs,
                               use_centroid = use_centroid)
    } else {

        pred_locs <- .SRE.predict(Sm=SRE_model,pred_locs=pred_locs,use_centroid=use_centroid)
    }
    pred_locs
}

setMethod("summary",signature(object="SRE"),
          function(object,...) {
              cat("SRE Object\n")
              cat("==========\n")
              cat("\n")
              cat(paste0("Formula: ",deparse(object@f)))
              cat("\n")
              cat(paste0("Number of datasets: ",length(object@data)))
              cat("\n")
              cat(paste0("Number of basis functions: ",object@basis@n))
              cat("\n")
              cat(paste0("Class of basis functions: ",class(object@basis)[1]))
              cat("\n")
              cat(paste0("Number of BAUs [extract using object@BAUs]: ",length(object@BAUs)))
              cat("\n")
              cat(paste0("Number of observations [extract using object@Z]: ",length(object@Z)))
              cat("\n")
              cat(paste0("Mean obs. variance at BAU level [extract using object@Ve]: ",mean(object@Ve@x)))
              cat("\n")
              cat(paste0("Fine-scale variatinal proportionality constant [extract using object@sigma2fshat]: ",object@sigma2fshat))
              cat("\n")
              cat(paste0("Dimensions of C in Z = C*Y + e [extract using object@Cmat]: ",deparse(dim(object@Cmat))))
              cat("\n")
              cat(paste0("Dimensions of S in Y = X*alpha + S*eta + delta [extract using object@S]: ",deparse(dim(object@S))))
              cat("\n")
              cat(paste0("Number of covariates: ",ncol(object@X)))
              cat("\n\n")
              cat(paste0("Summary of E(eta | Z) [extract using object@mu_eta]: \n"))
              cat("\n")
              print(summary(object@mu_eta[,1]))
              cat("\n\n")
              cat(paste0("Summary of Var(eta | Z) [extract using object@S_eta]: \n"))
              print(summary(diag(object@S_eta)))
              cat("\n\n")
              cat(paste0("Summary of Var(eta) [extract using object@Khat]: \n"))
              print(summary(diag(object@Khat)))
              cat("\n\n")
              cat(paste0("Regression coefficients [extract using object@alpha]: \n"))
              cat(deparse(as.vector(object@alphahat)))
})


.SRE.predict <- function(Sm,pred_locs,use_centroid) {
        depname <- all.vars(Sm@f)[1]
        pred_locs[[depname]] <- 0.1
        L <- .gstat.formula(Sm@f,data=pred_locs)
        X = as(L$X,"Matrix")
        if(is(pred_locs,"Spatial")) {
            if(use_centroid) {
                S0 <- eval_basis(Sm@basis,as.matrix(pred_locs[coordnames(Sm@data[[1]])]@data))
            } else {
                S0 <- eval_basis(Sm@basis,pred_locs)
            }
        } else if(is(pred_locs,"STFDF")) {
            if(use_centroid) {
                S0 <- eval_basis(Sm@basis,as.matrix(cbind(coordinates(pred_locs),pred_locs@data$t)))
            } else {
                stop("Can only use centroid when predicting with spatio-temporal data")
            }
        }
        pred_locs[[depname]] <- NULL

        alpha <- Sm@alphahat
        K <- Sm@Khat
        sigma2fs <- Sm@sigma2fshat
        D <- sigma2fs*Sm@Vfs + Sm@Ve
        Dinv <- Diagonal(x=1/diag(D))

        sig2_Vfs_pred <- Diagonal(x=sigma2fs*pred_locs$fs)

        if(is(pred_locs,"Spatial")) {
            idx <- match(row.names(pred_locs),row.names(Sm@BAUs))
        } else if (is(pred_locs,"STFDF")){
            idx <- match(pred_locs@data$n,Sm@BAUs@data$n)
        }
        C <- Sm@Cmat[,idx]

        if(sigma2fs >0) {
            LAMBDA <- as(bdiag(Sm@Khat,sig2_Vfs_pred),"symmetricMatrix")
            LAMBDAinv <- chol2inv(chol(LAMBDA))
            PI <- cBind(S0, .symDiagonal(n=length(pred_locs)))
            Qx <- t(PI) %*% t(C) %*% solve(Sm@Ve) %*% C %*% PI + LAMBDAinv
            temp <- cholPermute(Qx)
            ybar <- t(PI) %*%t(C) %*% solve(Sm@Ve) %*% (Sm@Z - C %*% X %*% alpha)
            x_mean <- cholsolve(Qx,ybar,perm=TRUE,cholQp = temp$Qpermchol, P = temp$P)
            Partial_Cov <- Takahashi_Davis(Qx,cholQp = temp$Qpermchol,P = temp$P)
            x_margvar <- diag(Partial_Cov)
            pred_locs[["mu"]] <- as.numeric(X %*% alpha + PI %*% x_mean)
        } else {
            LAMBDA <- as(Sm@Khat,"symmetricMatrix")
            LAMBDAinv <- chol2inv(chol(LAMBDA))
            PI <- S0
            Qx <- crossprod(solve(sqrt(Sm@Ve)) %*% C %*% PI) + LAMBDAinv
            #Qx <- t(PI) %*% t(C) %*% solve(Sm@Ve) %*% C %*% PI + LAMBDAinv
            ybar <- t(PI) %*%t(C) %*% solve(Sm@Ve) %*% (Sm@Z - C %*% X %*% alpha)
            Partial_Cov <- as(chol2inv(chol(Qx)),"dgeMatrix")  # Actually all Cov, convenient for later
            x_mean <- Partial_Cov %*% ybar
            x_margvar <- diag(Partial_Cov)
            pred_locs[["mu"]] <- as.numeric(X %*% alpha + PI %*% x_mean)
        }

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
            if(amp_factor > 1e12) {
                warning("sigma2fs is being estimated to zero. This might because because of an incorrect binnin procedure.")
                OK <- 1
            }
        }
        if(amp_factor > 1e12) {
            sigma2fs_new <- 0
            converged <- TRUE
        } else {
            sigma2fs_new <- uniroot(f = J,interval = c(sigma2fs/amp_factor,sigma2fs*amp_factor))$root
            D <- sigma2fs_new*Sm@Vfs + Sm@Ve
            Dinv <- Diagonal(x=1/diag(D))
            alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*% t(Sm@X) %*% Dinv %*% (Sm@Z - Sm@S %*% mu_eta)
            if(max(sigma2fs_new / sigma2fs, sigma2fs / sigma2fs_new) < 1.001) converged <- TRUE
        }
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
    if(is(BAUs,"Spatial"))
        if(!all(all.vars(f)[-1] %in% names(BAUs@data)))
            stop("All covariates need to be in the SpatialPolygons BAU object")
    if(is(BAUs,"ST"))
        if(!all(all.vars(f)[-1] %in% c(names(BAUs@data),coordnames(BAUs))))
            stop("All covariates need to be in the SpatialPolygons BAU object")
    if(!is(data,"list")) stop("Please supply a list of Spatial objects.")
    if(!all(sapply(data,function(x) is(x,"Spatial") | is(x,"ST")))) stop("All data list elements need to be of class Spatial or ST")
    if(!all(sapply(data,function(x) all.vars(f)[1] %in% names(x@data)))) stop("All data list elements to have values for the dependent variable")
    if(!all(sapply(data,function(x) identical(proj4string(x), proj4string(BAUs))))) stop("Please ensure all data items and BAUs have the same coordinate reference system")
    if(!(is(basis,"Basis") | is(basis,"TensorP_Basis"))) stop("basis needs to be of class Basis  or TensorP_Basis (package FRK)")
    #if(is(data,"SpatialPolygonsDataFrame") & !("std" %in% names(data))) stop("Polygon data needs to contain a field 'std' denoting the observation error")
    if(!("fs" %in% names(BAUs@data))) {
        warning("BAUs should contain a field 'fs' containing a basis function for fine-scale variation. Setting basis function equal to one everywhere.")
        BAUs$fs <- 1
    }
    if(!(is(BAUs,"SpatialPolygonsDataFrame") | is(BAUs,"STFDF"))) stop("BAUs should be a SpatialPolygonsDataFrame or a STFDF object")
    if(is(BAUs,"STFDF")) if(!is(BAUs@sp,"SpatialPolygonsDataFrame")) stop("The spatial component of the BAUs should be a SpatialPolygonsDataFrame")
    if((is(manifold(basis),"sphere")) & !all((coordnames(BAUs) == c("lon","lat")))) stop("Since a sphere is being used, please ensure that all coordinates (including those of BAUs) are in (lon,lat)")
    if(!est_error & !all(sapply(data,function(x) "std" %in% names(x@data)))) stop("If observational error is not going to be estimated, please supply a field 'std' in the data objects")
}

.gstat.formula <- function (formula, data)
{
    m = model.frame(terms(formula), as(data, "data.frame"), na.action = na.fail)
    Y = model.extract(m, "response")
    if (length(Y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)
    has.intercept = attr(Terms, "intercept")
    grid = numeric(0)
    xlevels = .getXlevels(Terms, m)
    list(y = Y, locations = coordinates(data), X = X, call = call,
         has.intercept = has.intercept, grid = as.double(unlist(grid)),
         xlevels = xlevels)
}
