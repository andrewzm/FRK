## Contains stuff for ICAR

#' @title Spatial Random Effects class
#' @description This is the central class definition of the \code{FRK} package, containing the model and all other information required for estimation and prediction.
#' @details The spatial random effects (SRE) model is the model employed in Fixed Rank Kriging, and the \code{SRE} object contains all information required for estimation and prediction from spatial data. Object slots contain both other objects (for example, an object of class \code{Basis}) and matrices derived from these objects (for example, the matrix \eqn{S}) in order to facilitate computations.
#'
#' @slot f formula used to define the SRE object. All covariates employed need to be specified in the object \code{BAUs}
#'@slot data the original data from which the model's parameters are estimated
#'@slot basis object of class \code{Basis} used to construct the matrix \eqn{S}
#'@slot BAUs object of class \code{SpatialPolygonsDataFrame}, \code{SpatialPixelsDataFrame} of \code{STFDF} that contains the Basic Areal Units (BAUs) that are used to both (i) project the data onto a common discretisation if they are point-referenced and (ii) provide a BAU-to-data relationship if the data has a spatial footprint
#' @slot S matrix constructed by evaluating the basis functions at all BAUs affected by the data (of class \code{Matrix})
#' @slot Ve measurement-error variance-covariance matrix (typically diagonal and of class \code{Matrix})
#' @slot Vfs fine-scale variance-covariance matrix at the data locations (typically diagonal and of class \code{Matrix}) up to a constant of proportionality estimated in the framework
#' @slot Vfs_BAUs fine-scale variance-covariance matrix at the BAU centroids (typically diagonal and of class \code{Matrix}) up to a constant of proportionality estimated in the framework
#' @slot Qfs_BAUs fine-scale precision matrix at the BAU centroids (typically diagonal and of class \code{Matrix}) up to a constant of proportionality estimated in the framework
#' @slot Z vector of observations (of class \code{Matrix})
#' @slot Cmat incidence matrix mapping the observations to the BAUs
#' @slot X matrix of covariates
#' @slot mu_eta updated expectation of random effects (estimated)
#' @slot S_eta updated covariance matrix of random effects (estimated)
#' @slot Khat prior covariance matrix of random effects (estimated)
#' @slot Khat_inv prior precision matrix of random effects (estimated)
#' @slot alphahat fixed-effect regression coefficients (estimated)
#' @slot sigma2fshat fine-scale variation scaling (estimated)
#' @keywords Spatial random effects, fixed rank kriging
setClass("SRE2",representation(data="list",
                              basis="Basis_obj",
                              BAUs="ANY",     # should be SpatialPolygonsDataFrame, SpatialPixelsDataFrame or STFDF
                              f = "formula",
                              S = "Matrix",
                              S0 = "Matrix",
                              Ve = "Matrix",
                              Vfs = "Matrix",
                              Vfs_BAUs = "Matrix",
                              Qfs_BAUs = "Matrix",
                              Z = "Matrix",
                              Cmat = "Matrix",
                              X = "Matrix",
                              mu_eta = "Matrix",
                              mu_xi = "Matrix",
                              S_eta = "Matrix",
                              Khat = "Matrix",
                              alphahat = "Matrix",
                              Q_eta = "Matrix",
                              Khat_inv = "Matrix",
                              B_run = "Matrix",
                              v_run = "Matrix",
                              sigma2fshat = "numeric",
                              fs_model = "character",
                              D_basis = "list",
                              K_type = "character",
                              lambda = "numeric"))

#' @title Construct SRE object, fit and predict
#' @description Main constructor of spatial random effects (SRE) object. Please see \code{\link{SRE-class}} for more details on the object's properties and methods.
#' @param f \code{R} formula relating the dependent variable (or transformations thereof) to covariates
#' @param data list of objects of class \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}
#' @param basis object of class \code{Basis}
#' @param BAUs object of class \code{SpatialPolygonsDataFrame}, the data frame which must contain covariate information as well as a field \code{fs} describing the fine-scale variation up to a constant of proportionality. If the function \code{FRK} is used directly, then BAUs are created automatically but only coordinates can then be used as covariates
#' @param est_error flag indicating whether the measurement-error variance should be estimated from variogram techniques. If this is set to 0, then \code{data} must contain a field \code{std}. Measurement-error estimation is not implemented for spatio-temporal datasets
#' @param average_in_BAU if \code{TRUE}, then multiple data points falling in the same BAU are averaged; the measurement error of the averaged data point is taken as the average of the individual measurement errors
#' @param fs_model if "ind" then the fine-scale variation is independent at the BAU level. If "ICAR", then an ICAR model is placed on the BAUs
#' @param SRE_model object returned from the constructor \code{SRE()}
#' @param n_EM maximum number of iterations for the EM algorithm
#' @param tol convergence tolerance for the EM algorithm
#' @param method parameter estimation method to employ. Currently only ``EM'' is supported
#' @param print_lik flag indicating whether likelihood should be printed or not on convergence of the estimation algorithm
# #' @param use_centroid flag indicating whether the basis functions are averaged over the BAU, or whether the basis functions are evaluated at the BAUs centroid in order to construct the matrix \eqn{S}. The flag can safely be set when the basis functions are approximately constant over the BAUs in order to reduce computational time
#' @param obs_fs flag indicating whether the fine-scale variation sits in the observation model (systematic error) or in the process model (process fine-scale variation, default)
#' @param pred_polys object of class \code{SpatialPoylgons} indicating the regions over which prediction will be carried out. The BAUs are used if this option is not specified
#' @param pred_time vector of time indices at which we wish to predict. All time points are used if this option is not specified
#' @param vgm_model an object of class \code{variogramModel} from the package \code{gstat} constructed using the function \code{vgm} containing the variogram model to fit to the data. The nugget is taken as the measurement error when \code{est_error = TRUE}. If unspecified the variogram used is \code{gstat::vgm(1, "Lin", d, 1)} where \code{d} is approximately one third of the maximum distance between any two points
#' @param K_type the parameterisation used for the \code{K} matrix. Currently this can be "unstructured" or "block-exponential"
#' @param lambda regularisation parameter (0 by default)
#' @param cross_validate the number \eqn{k} in \eqn{k}-fold cross-validation. If greater than 1, \code{lambda} is ignored and estimated through cross-validation
#' @param ... other parameters passed on to \code{auto_basis} and \code{auto_BAUs}
#' @details \code{SRE()} is the main function in the package as it constructs a spatial random effects model from the user-defined formula, data object, basis functions and a set of Basic Areal Units (BAUs). The function first takes each object in the list \code{data} and maps it to the BAUs -- this entails binning the point-referenced data into BAUs (and averaging within the BAU) if \code{average_in_BAU = TRUE}, and finding which BAUs are influenced by the polygon datasets. Following this, the incidence matrix \code{Cmat} is constructed, which appears in the observation model \eqn{Z = CY + e}, where \eqn{C} is the incidence matrix.
#'
#' The SRE model is given by \eqn{Y = T\alpha + S\eta + \delta}, where \eqn{T} are the covariates at BAU level, \eqn{\alpha} are the regression coefficients, \eqn{S} are the basis functions evaluated at the BAU level, \eqn{\eta} are the basis function weights, and \eqn{\delta} is the fine scale variation (at the BAU level). The covariance matrix of \eqn{\delta} is diagonal and proportional to the field `fs' in the BAUs (typically set to one). The constant of proportionality is estimated in the EM algorithm. All required matrices (\eqn{S,T} etc.) are computed and returned as part of the object, please see \code{\link{SRE-class}} for more details.
#'
#'\code{SRE.fit()} takes an object of class \code{SRE} and estimates all unknown parameters, namely the covariance matrix \eqn{K}, the fine scale variance \eqn{\sigma^2_{\delta}} and the regression parameters \eqn{\alpha}. The only method currently implemented is the Expectation Maximisation (EM) algorithm, which the user configures through \code{n_EM} and \code{tol}. The latter parameter, \code{tol}, is used as in Katzfuss and Cressie to, that is, the log-likelihood (given in Equation (16) in that work) is evaluated at each iteration at the current parameter estimate, and convergence is assumed to have been reached when this quantity stops changing by more than \code{tol}.
#'
#'The actual computations for the E-step and M-step are relatively straightforward. The E-step contains an inverse of an \eqn{n \times n} matrix, where \code{n} is the number of basis functions which should not exceed 2000. The M-step first updates the matrix \eqn{K}, which only depends on the sufficient statistics of the basis weights \eqn{\eta}. Then, the regression parameter \eqn{\alpha} is updated and a simple optimisation routine (a line search) is used to update the fine-scale variance \eqn{\sigma^2_{\delta}}. If the fine-scale errors and measurement errors are homoscedastic a closed-form solution is available for the update of \eqn{\sigma^2_{fs}}. Irrespectively, since the udpates of \eqn{\alpha} and \eqn{\sigma^2_{\delta}} are dependent, these two updates are iterated until the change in \eqn{\sigma^2_{\delta}} is no more than 0.1\%.
#'
#'Once the parameters are fitted, the \code{SRE} object is passed onto the function \code{SRE.predict()} in order to carry out optimal predictions over the same BAUs used to construct the SRE model with \code{SRE()}. The first part of the prediction process is to construct the matrix \eqn{S}. This is made computationally efficient by treating the prediction over polygons as that of the prediction over a combination of BAUs. This will yield valid results only if the BAUs are relatively small. Once the matrix \eqn{S} is found, a standard Gaussian inversion using the estimated parameters is used.
#'
#'\code{SRE.predict} returns the BAUs, which are of class \code{SpatialPolygonsDataFrame}, with two added attributes, \code{mu} and \code{var}. These can then be easily plotted using \code{spplot} or \code{ggplot2} (in conjunction with \code{\link{SpatialPolygonsDataFrame_to_df}}) as shown in the package vignettes.
#'\code{FRK}  runs \code{SRE}, \code{SRE.fit} and \code{SRE.predict} in successions with suitable defaults. It returns a list with the SRE object and the prediction polygons.
#' @references
#' Katzfuss, M., & Cressie, N. (2011). Spatio-temporal smoothing and EM estimation for massive remote-sensing data sets. Journal of Time Series Analysis, 32(4), 430--446.
#'
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
#' grid_BAUs <- auto_BAUs(manifold=real_line(),data=sim_data,
#'                        nonconvex_hull=FALSE,cellsize = c(0.01),type="grid")
#' grid_BAUs$fs = 1
#'
#' ### Set up SRE model
#' G <- auto_basis(manifold = real_line(),
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
#' S <- SRE.fit(S,n_EM = 5,tol = 0.01,print_lik=TRUE)
#'
#' ### Predict over BAUs
#' grid_BAUs <- SRE.predict(S)
#'
#' ### Plot
#' # X <- slot(grid_BAUs,"data") %>%
#' #      filter(x >= 0 & x <= 1)
#' # g1 <- LinePlotTheme() +
#' #    geom_line(data=X,aes(x,y=mu)) +
#' #    geom_errorbar(data=X,aes(x=x,ymax = mu + 2*sqrt(var), ymin= mu - 2*sqrt(var))) +
#' #    geom_point(data = data.frame(sim_data),aes(x=x,y=z),size=3) +
#' #    geom_line(data=sim_process,aes(x=x,y=proc),col="red")
#' # print(g1)
SRE2 <- function(f,data,basis,BAUs,est_error=FALSE,average_in_BAU = TRUE, fs_model = "ind",vgm_model = NULL, K_type = "block-exponential") {

    .check_args1(f=f,data=data,basis=basis,BAUs=BAUs,est_error=est_error)
    av_var <-all.vars(f)[1]
    ndata <- length(data)

    S <- Ve <- Vfs <- X <- Z <- Cmat <- list()

    print("Normalising basis function evaluations at BAU level ...")
    S0 <- eval_basis(basis,.polygons_to_points(BAUs))
    xx <- sqrt(rowSums((S0) * S0))
    xx <- xx + 1*(xx == 0) ## Where there are no basis functions do not divide zero by zero..
    S0 <- S0 / (xx %>% as.numeric())

    for(i in 1:ndata) {
        if(est_error) data[[i]]$std <- 0 ## Just set it to something, this will be overwritten later on
        if(est_error) {
            if(is(data[[i]],"ST"))
                stop("Estimation of error not yet implemented for spatio-temporal data")
            data_proc <- data[[i]]
            data_proc$Nobs <- 1
            data_proc <- est_obs_error(data_proc,variogram.formula=f, vgm_model = vgm_model)
            data[[i]]$std <- data_proc$std
        }


        print("Binning data ...")
        data_proc <- map_data_to_BAUs(data[[i]],
                                      BAUs,
                                      av_var = av_var,
                                      average_in_BAU = average_in_BAU)

        if(any(is.na(data_proc@data[av_var])))
            stop("NAs found when mapping data to BAUs. Do you have NAs in your data? If not, are you sure all your data are covered by BAUs?")

        L <- .extract.from.formula(f,data=data_proc)
        X[[i]] <- as(L$X,"Matrix")
        Z[[i]] <- Matrix(L$y)
        Ve[[i]] <- Diagonal(x=data_proc$std^2)


        C_idx <- BuildC(data_proc,BAUs)

        Cmat[[i]] <- sparseMatrix(i=C_idx$i_idx,
                                  j=C_idx$j_idx,
                                  x=1,
                                  dims=c(length(data_proc),
                                         length(BAUs)))

        if(any(rowSums(Cmat[[i]])==0))
            stop("I have found difficulty in associating the data with the BAUs.
                 If you have point-referenced data
                 then this could be because you have data outside BAUs. If you have
                 polygon data, then this could be because no BAUs centroids are
                 within the polygons. For polygon data, influence on a BAU is determined from
                 whether the BAU centroid falls within the polygon or not.")

        Cmat[[i]] <- Cmat[[i]] / rowSums(Cmat[[i]]) ## Average BAUs for polygon observations

        if(fs_model == "ind") {
            Vfs[[i]] <- tcrossprod(Cmat[[i]] %*% Diagonal(x=sqrt(BAUs$fs)))
        } else if(fs_model == "ICAR") {
            Vfs[[i]] <- Matrix(ncol=0,nrow=0) ## Ignore variance matrix if ICAR
        } else {
            stop("Model needs to be ``ind'' or ``ICAR''.")
        }


        # if(length(Cmat[[i]]@x) == nrow(Cmat[[i]]))  { # if point observations
        #     Vfs[[i]] <- Diagonal(x = Cmat[[i]] %*% sqrt(BAUs$fs) %>% as.numeric())
        # } else {
        #     Vfs[[i]] <- Diagonal(x = rep(0,nrow(Cmat[[i]])))
        # }

        ## The following code was used when we assumed the BAUs were large
        ## compared to basis functions, we now have deprecated it
        # print("Evaluating basis functions at observation locations...")
        # S[[i]] <- eval_basis(basis, s = data_proc)
        # print("Done.")

        S[[i]] <- Cmat[[i]] %*% S0

        ## Note that S constructed in this way is similar to Cmat %*% S_BAUs where S_BAUs is the
        ## basis functions evaluated at the BAUs. Verify this by checking the following are similar
        ## S2 <- eval_basis(basis, s = BAUs)
        ## (Cmat[[i]] %*% S2) - S[[i]]
    }

    if(fs_model == "ind") {
        Qfs_BAUs <- Diagonal(x=1/BAUs$fs)
        Vfs_BAUs <- Diagonal(x=BAUs$fs)
    } else if(fs_model == "ICAR") {
        ## Make block diagonal for spatio-temporal
        message("Finding the polygon neighbours...")
        ## Caters for both spatial and ST
        nblist <- spdep::poly2nb(as(BAUs,"SpatialPolygonsDataFrame")[,1][,1])
        Qfs_BAUs <- .prec_from_neighb(nblist)
        if(is(BAUs,"STFDF")) {
            Qfs_BAUs <- do.call("bdiag",
                                lapply(1:length(BAUs@time),function(i) {Qfs_BAUs}))
        }
        Vfs_BAUs <- Matrix(ncol=0,nrow=0)
    }

    S <- do.call("rBind",S)
    X <- do.call("rBind",X)
    Z <- do.call("rBind",Z)
    Ve <- do.call("bdiag",Ve)
    Vfs <- do.call("bdiag",Vfs)
    D_basis <- BuildD(basis)
    #K_norm <- .initialise_K(basis,D_basis)
    #K_init <- var(Z[,1])*K_norm
    K_init = Diagonal(n=nbasis(basis),x = 1/(1/var(Z[,1])))
    K_inv_init = Diagonal(n=nbasis(basis),x = (1/var(Z[,1])))

    new("SRE",
        data=data,
        basis=basis,
        BAUs=BAUs,
        f = f,
        S = S,
        S0 = S0,
        Ve = Ve,
        Vfs = Vfs,
        Vfs_BAUs = Vfs_BAUs,
        Qfs_BAUs = Qfs_BAUs,
        Z = Z,
        Cmat = do.call("rBind",Cmat),
        X = X,
        mu_eta = Matrix(0,nbasis(basis),1),
        S_eta = Diagonal(x = rep(1,nbasis(basis))),
        Q_eta = Diagonal(x = rep(1,nbasis(basis))),
        Khat = K_init,
        Khat_inv = K_inv_init,
        alphahat = solve(t(X) %*% X) %*% t(X) %*% Z,
        sigma2fshat = mean(diag(Ve)) / 4,
        B_run = Diagonal(n=nbasis(basis),x = 1/var(Z[,1])),
        v_run = Matrix(0,nbasis(basis),nbasis(basis)),
        fs_model = fs_model,
        D_basis = D_basis,
        K_type = K_type,
        lambda = 0)
}

.SRE.EMstep.ICAR <- function(Sm) {

    alpha <- Sm@alphahat
    K <- Sm@Khat
    Kinv <- Sm@Khat_inv
    sigma2fs <- Sm@sigma2fshat
    Qfs_norm <- Sm@Qfs_BAUs %>% as("dgTMatrix")
    Cmat <- Sm@Cmat
    r <- nrow(K)
    n <- length(Sm@BAUs)
    Qe <- solve(Sm@Ve)

    GAMMA <- as(bdiag(Kinv,(1/sigma2fs) * Qfs_norm),"symmetricMatrix") %>% as("dgTMatrix")
    PI <- cBind(Sm@S, Cmat %*% .symDiagonal(n=length(Sm@BAUs)))
    Qx <- (t(PI) %*% solve(Sm@Ve) %*% PI + GAMMA) %>% as("dgTMatrix")

    ## Add (zero) elements to Qx so that all covariance elements associated with eta are computed
    ## This may be removed in the future if we work with uniformly sparse K
    ij <- expand.grid(i=0:(r-1),j=0:(r-1))
    miss_idx <- setdiff(ij,data.frame(i=Qx@i,j=Qx@j))
    Qx@i <- c(Qx@i,miss_idx[,1])
    Qx@j <- c(Qx@j,miss_idx[,2])
    Qx@x <- c(Qx@x,rep(0L,nrow(miss_idx)))
    Qx <- as(Qx,"dgCMatrix")
    temp <- cholPermute(Qx)
    ybar <- t(PI) %*% Qe %*% (Sm@Z - Sm@X %*% alpha)
    x_mean <- cholsolve(Qx,ybar,perm=TRUE,cholQp = temp$Qpermchol, P = temp$P)
    Cov <- Takahashi_Davis(Qx,cholQp = temp$Qpermchol,P = temp$P) # PARTIAL

    MeanOuter_sparse <- sparseMatrix(i=GAMMA@i + 1, j=GAMMA@j + 1,
                                     x = x_mean[GAMMA@i+1] * x_mean[GAMMA@j+1])
    MeanOuter_eta <- tcrossprod(Matrix(x_mean[1:r]))

    Sm@Khat <- .regularise_K(Sm = Sm,
                             S_eta = as(forceSymmetric(Cov[(1:r),(1:r)]),"symmetricMatrix"),
                             mu_eta = (Matrix(x_mean[1:r])))
    Sm@Khat_inv <- chol2inv(chol(Sm@Khat))

    Sm@sigma2fshat <- sum(Qfs_norm * (Cov[-(1:r),-(1:r)] + MeanOuter_sparse[-(1:r),-(1:r)]))/ length(Sm@BAUs)
    Sm@alphahat <- solve(t(Sm@X) %*% Qe %*% Sm@X) %*% t(Sm@X) %*% Qe %*% (Sm@Z - PI %*% x_mean)
    Sm@mu_eta <- Matrix(x_mean[1:r])
    Sm@mu_xi <- Matrix(x_mean[-(1:r)])
    Sm@S_eta <- Cov[1:r,1:r]
    Sm

}

.loglik.ICAR <- function(Sm) {

    # warning("Monitoring complete-data likelihood")
    # res <- Sm@Z - Sm@X %*% Sm@alphahat - Sm@S %*% Sm@mu_eta - Sm@Cmat %*% Sm@mu_xi
    # (-0.5 * t(res) %*% solve(Sm@Ve) %*%  res) %>% as.numeric()
    S <- Sm@S
    K <- Sm@Khat
    chol_K <- chol(K)
    Kinv <- chol2inv(chol_K)
    resid <- Sm@Z - Sm@X %*% Sm@alphahat
    N <- length(Sm@Z)
    Qe <- solve(Sm@Ve)
    Cmat <- Sm@Cmat
    Qfs <- (1/Sm@sigma2fshat) * Sm@Qfs_BAUs
    R <- chol(Qfs  + t(Cmat) %*% Qe %*% Cmat)
    Dinv <- Qe*1.000000001 - crossprod(t(solve(R)) %*% t(Cmat) %*% Qe)
    chol_Dinv <- chol(Dinv)
    D <- chol2inv(chol_Dinv)
    cholD <- chol(D)
    cholDinvT <- t(solve(cholD))
    S_Dinv_S <-  crossprod(cholDinvT %*% S)
    log_det_SigmaZ <- determinant(Kinv + S_Dinv_S,logarithm = TRUE)$modulus +
        determinant(K,logarithm = TRUE)$modulus +
        logdet(cholD)

    ## Alternatively: (slower but more direct)
    # Dinv <- chol2inv(chol(D))
    # SigmaZ_inv <- Dinv - Dinv %*% S %*% solve(Kinv + S_Dinv_S) %*% t(S) %*% Dinv
    # SigmaZ_inv2 <- Dinv - tcrossprod(Dinv %*% S %*% solve(R))

    R <- chol(Kinv + S_Dinv_S)

    rDinv <- crossprod(cholDinvT %*% resid,cholDinvT)
    ## Alternatively: # rDinv <- t(resid) %*% Dinv

    quad_bit <- crossprod(cholDinvT %*% resid) - tcrossprod(rDinv %*% S %*% solve(R))
    ## Alternatively: # quad_bit <- rDinv %*% resid - tcrossprod(rDinv %*% S %*% solve(R))

    llik <- -0.5 * N * log(2*pi) -
        0.5 * log_det_SigmaZ -
        0.5 * quad_bit
    as.numeric(llik)
}


#### Using cross-validation parameter

#' @rdname SRE
#' @export
SRE.fit2 <- function(SRE_model,n_EM = 100L, tol = 0.01, lambda = 0, method="EM", print_lik=FALSE, cross_validate=1L) {
    .check_args2(n_EM = n_EM,tol = tol,lambda = lambda,method = method,print_lik = print_lik,cross_validate = cross_validate)
    if(!(length(cross_validate) == 1 | length(cross_validate) == nrow(count_res(SRE_model))))
        stop("cross_validate needs to be of length one or of length equal to the number of basis-function resolutions")

    if(any(cross_validate > 1L) & length(lambda) == 1) {
        VarZ <- var(SRE_model@Z[,1])
        lambda = c(0,1/VarZ,10/VarZ,100/VarZ)
        print("Cross-validating on lambda = ")
        print(lambda)
    }


    if(all(cross_validate == 1)) {
        .SRE.fit(SRE_model = SRE_model, n_EM = n_EM, tol = tol, lambda = lambda, method = "EM", print_lik=print_lik)
    } else {
        m <- nrow(SRE_model@Z)
        all_coords <- SRE_model@Cmat %*% coordinates(SRE_model@BAUs)

        nres <- nrow(count_res(SRE_model))
        current_lambda <- rep(0,nres)
        max_lambda = Inf  ## Make sure lambdas are monotonic in resolution
        if(length(cross_validate) == 1) {
            nres = 1 # Only optimise one lambda (for all resolutions)
        }
        for(res in nres : 1) {
            num_at_res <- count_res(SRE_model)[res,]$n
            if(num_at_res < m/4 & length(cross_validate) > 1) {
                print(paste0("Dividing data into ",num_at_res*2," clusters"))
                numclusters = Hmisc::ceil(min(num_at_res*2, m))
                cluster_labels <- kmeans(all_coords,centers = numclusters)$cluster
            } else {
                numclusters = m
                cluster_labels <- sample(1:m,m)
            }

            partitions <- cut(1:numclusters,
                              seq(0,numclusters+1,length=cross_validate[nres] + 1),
                              labels=FALSE)
            ESS <- crps <- sq_resid <- ESS_score <- crps_score <- cv_score <- NULL

            for(l in seq_along(lambda[lambda <= max_lambda])) {
                if(length(cross_validate) > 1) {
                    current_lambda[res] <- lambda[l]
                } else {
                    current_lambda <- rep(lambda[l],nres)
                }
                for (i in 1:cross_validate[nres]) {
                    these_clusters <- which(partitions == i)
                    rm_idx <- which(cluster_labels %in% these_clusters)
                    S_part <- .remove_obs_from_SRE(S = SRE_model, rm_idx = rm_idx)
                    S_part <- .SRE.fit(SRE_model = S_part, n_EM = n_EM, tol = tol, lambda = current_lambda,
                                       method = "EM", print_lik=FALSE)
                    BAUs_to_predict <- apply(SRE_model@Cmat[rm_idx,],1,function(x) which(x==1))
                    Validate_obs <- SRE.predict(SRE_model = S_part,      # SRE model
                                                pred_polys = S_part@BAUs[BAUs_to_predict,],
                                                obs_fs = FALSE)
                    sq_resid[i] <- mean((Validate_obs$mu - SRE_model@Z[rm_idx])^2)
                    crps[i] <- verification::crps(SRE_model@Z[rm_idx],
                                                  cbind(Validate_obs$mu,sqrt(Validate_obs$var + diag(SRE_model@Ve)[rm_idx])))$CRPS
                    ESS[i] <- mean((Validate_obs$var + diag(SRE_model@Ve)[rm_idx] -
                                        (Validate_obs$mu - SRE_model@Z[rm_idx])^2)^2)
                    #hist((Validate_obs$mu - SRE_model@Z[rm_idx])/sqrt(Validate_obs$var + diag(SRE_model@Ve)[rm_idx]))
                }
                cv_score[l] <- mean(sq_resid)
                crps_score[l] <- mean(crps)
                ESS_score[l] <- mean(ESS)
            }
            print(paste0("Cross validation results for ",
                         ifelse(length(cross_validate) > 1,paste0("resolution ",res),"all resolutions"),":"))
            print("---------------------------------------------")
            print(data.frame(lambda = lambda[lambda <= max_lambda],
                             sq_res = cv_score, crps=crps_score, ESS = ESS_score))
            lambda_best <- lambda[which.min(cv_score)]
            print(paste0("Proceeding with lambda = ", lambda_best," for this(these) resolution(s)"))
            current_lambda[res] <- lambda_best
            max_lambda <- lambda_best
        }

        .SRE.fit(SRE_model = SRE_model, n_EM = n_EM, tol = tol, lambda = current_lambda, method = "EM", print_lik=print_lik)

    }

}

.check_args2B <- function(n_EM = 100L, tol = 0.01, lambda = 0, method="EM", print_lik=FALSE,...) {
    if(!is.numeric(n_EM)) stop("n_EM needs to be an integer")
    if(!(n_EM <- round(n_EM)) > 0) stop("n_EM needs to be greater than 0")
    if(!is.numeric(tol)) stop("tol needs to be a number greater than zero")
    if(!(tol > 0)) stop("tol needs to be a number greater than zero")
    if(!(method == "EM")) stop("Currently only the EM algorithm is implemented for parameter estimation")
    if(!(is.logical(print_lik))) stop("print_lik needs to be a logical quantity")
    if(!(is.numeric(lambda))) stop("lambda needs to be a number")
    if(!(all(lambda >= 0))) stop("lambda needs to be greater or equal to zero")
    if(!(is.integer(cross_validate))) stop("cross_validate needs to be an integer or a vector of integers")
    if(length(lambda) > 1 & cross_validate == 1L) stop("to find an optimal lambda cross_validate must be greater than one to split the data into training sets and validation sets")
}


.SRE.fit2 <- function(SRE_model,n_EM = 100L, tol = 0.01, lambda = 0, method="EM", print_lik=FALSE) {
    n <- nbasis(SRE_model)
    X <- SRE_model@X
    lk <- rep(0,n_EM)
    SRE_model@lambda <- lambda

    if(!is.na(tol) & (SRE_model@fs_model == "ICAR")) {
        warning("Cannot monitor the observed likelihood with the ICAR model.
                Monitoring changes in eta instead.")
        lik_plot_ylab <- "norm(eta)"
    } else {
        lik_plot_ylab <- "log likelihood"
    }

    if(opts_FRK$get("progress")) pb <- utils::txtProgressBar(min = 0, max = n_EM, style = 3)
    for(i in 1:n_EM) {
        if (!(SRE_model@fs_model == "ICAR")){
            #print(system.time( lk[i] <- .loglik(SRE_model)))  # Compute likelihood
            lk[i] <- .loglik(SRE_model)
        } else {
            lk[i] <- sqrt(sum(SRE_model@mu_eta^2))
        }

        # Still in development:
        # if(i == 2) {
        #     print("Normalising basis function evaluations at BAUs ...")
        #     S0 <- eval_basis(SRE_model@basis,as.matrix(SRE_model@BAUs[coordnames(SRE_model@data[[1]])]@data))
        #     xx <<- sqrt(rowSums((S0 %*% SRE_model@S_eta) * S0))
        #     S0 <- S0/xx
        #     SRE_model@S <- SRE_model@Cmat %*% S0
        #     print("Done ...")
        # }

        SRE_model <- .SRE.Estep(SRE_model)
        SRE_model <- .SRE.Mstep(SRE_model)
        if(opts_FRK$get("progress")) utils::setTxtProgressBar(pb, i)
        if(i>1)
            if(abs(lk[i] - lk[i-1]) < tol) {
                print("Minimum tolerance reached")
                break
            }
    }
    if(opts_FRK$get("progress")) close(pb)
    #if(SRE_model@sigma2fshat == 0)
    #warning("sigma2fs is being estimated to zero.
    #This might because of an incorrect binning procedure or because
    #too much measurement error is being assumed (or because the latent
    #field is indeed that smooth, but unlikely).")

    if(i == n_EM) print("Maximum EM iterations reached")
    if(print_lik & !is.na(tol)) {
        plot(1:i,lk[1:i],ylab=lik_plot_ylab,xlab="EM iteration")
    }
    SRE_model
}

.SRE.Estep2 <- function(Sm) {
    if(Sm@fs_model == "ind") {
        Sm <- .SRE.Estep.ind(Sm)
    } else if(Sm@fs_model == "ICAR") {
        Sm <- .SRE.EMstep.ICAR(Sm)
    }
}

.SRE.Mstep2 <- function(Sm) {
    if(Sm@fs_model == "ind") {
        Sm <- .SRE.Mstep.ind(Sm)
    } else if(Sm@fs_model == "ICAR") {
        Sm # M-step already carried out
    }
}

.SRE.Estep.ind2 <- function(Sm) {
    alpha <- Sm@alphahat
    K <- Sm@Khat
    Kinv <- Sm@Khat_inv
    sigma2fs <- Sm@sigma2fshat

    D <- sigma2fs*Sm@Vfs + Sm@Ve
    if(is(D,"dtCMatrix")) {
        cholD <- sqrt(D)
        cholDinv <- solve(cholD)
        Dinv <- solve(D)
    } else {
        cholD <- Matrix::chol(D)
        cholDinv <- solve(cholD)
        Dinv <- chol2inv(chol(D))
    }

    Q_eta <- (crossprod(t(cholDinv) %*% Sm@S) + Kinv)
    S_eta <- chol2inv(chol(Q_eta))
    mu_eta <- (S_eta) %*%(t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha))


    ## Deprecated:
    # if(!is(Q_eta,"dsCMatrix")) Q_eta <- as(Q_eta,"dsCMatrix")
    # chol_Q_eta <- cholPermute(Q_eta)
    # mu_eta <- cholsolve(Q_eta,(t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha)),
    #                     perm=TRUE, cholQp = chol_Q_eta$Qpermchol,P = chol_Q_eta$P)
    # S_eta <- Matrix()

    ## Deprecated:
    # S_eta <- chol2inv(chol(crossprod(t(cholDinv) %*% Sm@S) + Kinv))
    # mu_eta <- S_eta %*% (t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha))

    Sm@mu_eta <- mu_eta
    Sm@S_eta <- S_eta
    Sm@Q_eta <- Q_eta
    Sm

}

########################################
####### DEPRECATED #####################
########################################

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
#' \dontrun{library(ggplot2)
#'    show_basis(G_space_time)}
#' @export
sp_to_ST_basis <- function(G_spatial,t_knots = 1,manifold=STsphere()) {

    if(!dimensions(manifold(G_spatial))==2)
        stop("Please ensure that the dimension of the spatial manifold is 2")
    if(!dimensions(manifold)==3)
        stop("Please ensure that the dimension of the space-time manifold is 3")
    if(!is.numeric(t_knots))
        stop("Please ensure the knots are of type numeric")

    n <- nbasis(G_spatial)  # extract the number of basis functions
    G <- list()             # initialise the new basis list (one list item per time point)
    for(i in seq_along(t_knots)) { # for each time point
        Gt <- G_spatial            # create a new set of basis functions identical to the spatial ones
        sapply(1:n, function(j) {  # for each basis function
            this_c <-   get("c",environment(Gt@fn[[j]]))    # retrieve centroid
            new_c <- cbind(this_c,t_knots[i])               # add time coordinate to centroid
            assign("c",new_c,environment(Gt@fn[[j]]))       # assign new coordinate to function environment
        })
        Gt@df <- cbind(Gt@df,loc3=t_knots[i])   # update the data-frame with the new time coordinate
        Gt@manifold <- manifold                 # put basis function on the space-time manifold
        G[[i]] <- Gt                            # assign to this time point
    }
    G <- Reduce("concat",G)  # concatenate using the S4 function concat
}

## Evaluate basis over BAUs... deprecated?
.eval_basis.BAUs <- function(basis,BAUs,use_centroid) {
    if(is(BAUs,"Spatial")) {
        if(use_centroid) {
            #S0 <- eval_basis(Sm@basis,as.matrix(BAUs[coordnames(Sm@data[[1]])]@data))
            S0 <- eval_basis(basis,.polygons_to_points(BAUs))
        } else {
            S0 <- eval_basis(basis,BAUs)
        }
    } else if(is(BAUs,"STFDF")) {
        if(use_centroid) {
            #S0 <- eval_basis(Sm@basis,as.matrix(cbind(coordinates(BAUs),BAUs@data$t)))
            S0 <- eval_basis(basis,.polygons_to_points(BAUs))
        } else {
            stop("Can only use centroid when predicting with spatio-temporal data")
        }
    }
    S0
}


#' @name plotting-themes
#' @aliases LinePlotTheme
#' @aliases EmptyTheme
#' @title Plotting themes
#' @description Formats a ggplot object for neat plotting.
#' @return Object of class \code{ggplot}
#' @keywords ggplot
#' @export
#' @details \code{LinePlotTheme()} creates \code{ggplot} object with a white background, a relatively large font and grid lines. \code{EmptyTheme()} on the other hand creates a \code{ggplot} object with no axes or legends.
#' @examples
#' \dontrun{
#' X <- data.frame(x=runif(100),y = runif(100), z = runif(100))
#' LinePlotTheme() + geom_point(data=X,aes(x,y,colour=z))
#' EmptyTheme() + geom_point(data=X,aes(x,y,colour=z))}

#' @rdname plotting-themes
#' @export
LinePlotTheme <- function() {
    g <- ggplot() + theme(panel.background = element_rect(fill='white', colour='black'),text = element_text(size=20),
                          panel.grid.major =  element_line(colour = "light gray", size = 0.05),
                          panel.border  = element_rect(fill=NA, colour='black'))
    #plot.margin=unit(c(5,5,5,0),"mm"))
    return (g)
}

#' @rdname plotting-themes
#' @export
EmptyTheme <- function() {
    g <- ggplot() +  theme(panel.background = element_rect(fill='white', colour='white'),panel.grid=element_blank(),axis.ticks=element_blank(),
                           panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
    return (g)
}


clip_polygons_lonlat <- function(d,key) {
    lon <- lat <- NULL
    plyr::ddply(d,key,function(df) {
        if(diff(range(df$lon)) > 90) {
            Y1 <- filter(df,lon >= 0) %>% mutate(id= df[key][1,]*1e6)
            Y1$lon[which(Y1$lon %in% sort(Y1$lon,decreasing=T)[1:2])] <- 179.99
            Y1 <- rbind(Y1,Y1[1,])
            Y2 <- filter(df,lon <= 0) %>% mutate(id= df[key][1,]*1e6+1)
            Y2$lon[which(Y2$lon %in% sort(Y2$lon,decreasing=F)[1:2])] <- -179.99
            Y2 <- rbind(Y2,Y2[1,])
            rbind(Y1,Y2)
        } else {
            df
        }})
}


.initialise_K <- function(basis,D_basis) {

    all_res <- count_res(basis)
    K_norm <- lapply(1:nrow(all_res),
                     function(i) {
                         idx <- which(basis@df$res == i)
                         if(is(basis,"TensorP_Basis")) {
                             tau_init_space <- max(D_basis$Basis1[[i]])/5 ## space
                             tau_init_time <- 3 ## time
                             Ki <- kronecker(exp(-D_basis$Basis2[[1]]/tau_init_time),
                                             exp(-D_basis$Basis1[[i]]/tau_init_space))
                         } else {
                             tau_init_space <- max(D_basis[[i]])/5 ## space
                             Ki <- exp(-D_basis[[i]]/tau_init_space)
                         }
                     })

    K_norm <- do.call("bdiag",K_norm)
    idx_all <- unlist(lapply(1:nrow(all_res), function(i) which(basis@df$res == i)))

    # Rearrange in order time/resolution when we have tensor products
    # When we don't have tensor product idx_all and 1:nrow(K) should be the same
    K_norm <- reverse_permute(K_norm,idx_all)


}

.remove_obs_from_SRE <- function(S, rm_idx) {
    S_part <- S
    S_part@S <- S@S[-rm_idx,]
    S_part@Ve <- S@Ve[-rm_idx,-rm_idx]
    S_part@Vfs <- S@Vfs[-rm_idx,-rm_idx]
    S_part@Z <- S@Z[-rm_idx,,drop=FALSE]
    S_part@Cmat <- S@Cmat[-rm_idx,,drop=FALSE]
    S_part@X <- S@X[-rm_idx,,drop=FALSE]
    S_part
}


## The following function is the internal prediction function
.SRE.predict2 <- function(Sm,obs_fs = FALSE,pred_polys = NULL,pred_time = NULL) {

    ## If the user does not specify time points to predict at when in space-time
    ## Then predict at every time point
    if(is.null(pred_time) & is(Sm@BAUs,"ST"))
        pred_time <- 1:length(Sm@BAUs@time)

    ## We start by assuming that we will predict at BAUs
    predict_BAUs <- TRUE

    ## Get BAUs from the SRE model
    BAUs <- Sm@BAUs

    ## If the user has not specified polygons over which to predict, then CP is
    ## just the diagonal matrix and we predict over all the BAUs
    if(is.null(pred_polys)) {
        CP <- Diagonal(length(BAUs))
    } else {
        ## The user has maybe specified a subset of (could be all) the BAUs over which to predict.
        ## The following checks whether pred_polys is a subset of the BAUs through the row names
        pred_polys_are_BAUs <- all(row.names(pred_polys) %in% row.names(BAUs))

        ## If the user has specified a subset of BAUs
        if(pred_polys_are_BAUs) {
            ## See which BAUs the user has specified
            BAUs_idx <- match(row.names(pred_polys), row.names(BAUs))

            ## Construct an incidence matrix that picks out these BAUs
            CP <-  sparseMatrix(i=1:length(pred_polys),
                                j=BAUs_idx,
                                x=1,
                                dims=c(length(pred_polys),
                                       length(BAUs)))

        } else {
            ## The user has specified arbitrary polygons
            ## First try to coerce what the user supplied to Polygons (not pixels etc.)
            ## Recall that for now only Spatial pred_polys are allowed so the following is
            ## always valid
            pred_polys <- as(pred_polys,"SpatialPolygonsDataFrame")

            ## Based on these polygons construct the C matrix
            C_idx <- BuildC(pred_polys,BAUs)
            CP <- sparseMatrix(i=C_idx$i_idx,
                               j=C_idx$j_idx,
                               x=1,
                               dims=c(length(pred_polys),
                                      length(BAUs)))

            ## As in SRE(), make sure the polgons are averages (not sums)
            CP <- CP / rowSums(CP)

            ## If even one polygon encompasses more than one BAU, then we need to
            ## predict over linear combinations of BAUs, and hence need to
            ## compute the full covariance matrix. Note this by setting
            ## predict_BAUs <- FALSE
            if(!all(table(C_idx$i_idx) == 1))
                predict_BAUs <- FALSE   ## Need to compute full covariance matrix
        }
    }

    ## Get the CZ matrix
    CZ <- Sm@Cmat

    ## If the user has specified which polygons he want we can remove the ones we don't need
    ## We only need those BAUs that are influenced by observations and prediction locations
    if(!is.null(pred_polys)) {

        ## The needed BAUs are the nonzero column indices of CZ and CP
        needed_BAUs <- union(as(CP,"dgTMatrix")@j+1, as(CZ,"dgTMatrix")@j+1)

        ## Filter the BAUs and the matrices
        BAUs <- BAUs[needed_BAUs,]
        CP <- CP[,needed_BAUs]
        CZ <- CZ[,needed_BAUs]
        Sm@S0 <- Sm@S0[needed_BAUs,]
    }

    # Deprecated:
    # if(is(BAUs,"ST")){
    #     needed_BAUs <- BAUs[,pred_time]$n
    #     BAUs <- BAUs[,pred_time]
    #     CP <- CP[,needed_BAUs]
    #     CZ <- CZ[,needed_BAUs]
    # }

    ## Retrieve the dependent variable name
    depname <- all.vars(Sm@f)[1]

    ## Set the dependent variable in BAUs to something just so that .extract.from.formula doesn't
    ## throw an error.. we will NULL it shortly after
    BAUs[[depname]] <- 0.1

    ## Extract covariates from BAUs
    L <- .extract.from.formula(Sm@f,data=BAUs)
    X = as(L$X,"Matrix")
    BAUs[[depname]] <- NULL

    ## Set variables to make code more concise
    S0 <- Sm@S0
    alpha <- Sm@alphahat
    K <- Sm@Khat
    sigma2fs <- Sm@sigma2fshat
    mu_eta <- Sm@mu_eta
    S_eta <- Sm@S_eta

    if(Sm@fs_model == "ind") {
        D <- sigma2fs*Sm@Vfs + Sm@Ve
        if(is(D,"dtCMatrix")) {
            Dchol <- sqrt(D)
            Dinv <- solve(D)
        } else {
            Dchol <- chol(D)
            Dinv <- chol2inv(Dchol)
        }

        sig2_Vfs_pred <- Diagonal(x=sigma2fs*BAUs$fs)
        Q <- solve(sig2_Vfs_pred)


    } else if(Sm@fs_model == "ICAR") {
        Q <- (1/sigma2fs) * Sm@Qfs_BAUs
    }

    if(is(BAUs,"Spatial")) {
        idx <- match(row.names(BAUs),row.names(Sm@BAUs))
    } else if (is(BAUs,"STFDF")){
        idx <- match(BAUs@data$n,Sm@BAUs@data$n)
    }

    if(!obs_fs) {
        if(sigma2fs >0) {
            #LAMBDA <- as(bdiag(Sm@Khat,sig2_Vfs_pred),"symmetricMatrix")
            LAMBDAinv <- bdiag(Sm@Khat_inv,Q)
            PI <- cBind(S0, .symDiagonal(n=length(BAUs)))
            tC_Ve_C <- t(CZ) %*% solve(Sm@Ve) %*% CZ + 0*.symDiagonal(ncol(CZ)) # Ensure zeros
            Qx <- t(PI) %*% tC_Ve_C %*% PI + LAMBDAinv
            chol_Qx <- cholPermute(as(Qx,"dgCMatrix"))
            ybar <- t(PI) %*%t(CZ) %*% solve(Sm@Ve) %*% (Sm@Z - CZ %*% X %*% alpha)
            x_mean <- cholsolve(Qx,ybar,perm=TRUE,cholQp = chol_Qx$Qpermchol, P = chol_Qx$P)
            if(predict_BAUs) {
                Cov <- Takahashi_Davis(Qx,cholQp = chol_Qx$Qpermchol,P = chol_Qx$P) # PARTIAL
                BAUs[["var"]] <- .batch_compute_var(S0,Cov,obs_fs = !(!obs_fs & sigma2fs > 0))
                BAUs[["sd"]] <- sqrt(BAUs[["var"]])
            } else {
                ## Do not compute covariance now
                #Cov <- cholsolve(Qx,Diagonal(nrow(Qx)),perm=TRUE,
                #                 cholQp = chol_Qx$Qpermchol, P = chol_Qx$P) # FULL
            }


        } else {
            LAMBDA <- as(Sm@Khat,"symmetricMatrix")
            LAMBDAinv <- chol2inv(chol(LAMBDA))
            PI <- S0
            Qx <- crossprod(solve(sqrt(Sm@Ve)) %*% CZ %*% PI) + LAMBDAinv
            ybar <- t(PI) %*%t(CZ) %*% solve(Sm@Ve) %*% (Sm@Z - CZ %*% X %*% alpha)
            Cov <- as(chol2inv(chol(Qx)),"dgeMatrix")  ## Do all covariance matrix
            ## We can do all the covariance matrix since the dimension is equal to those of eta
            x_mean <- Cov %*% ybar
            ## variance too hard to compute all at once -- do it in blocks of 1000
            BAUs[["var"]] <- .batch_compute_var(S0,Cov,obs_fs = !(!obs_fs & sigma2fs > 0))
            BAUs[["sd"]] <- sqrt(BAUs[["var"]])
        }
        BAUs[["mu"]] <- as.numeric(X %*% alpha + PI %*% x_mean)

        ### Since we have all the elements we can use first principles from the sparse covariance matrix
        #BAUs[["var"]] <- .batch_compute_var(PI,Cov)

    }

    if(obs_fs) {
        Qobs <- solve(Sm@Ve)
        Qx <- (crossprod(t(sqrt(Qobs)) %*% (Sm@S %>% as("dgCMatrix"))) + chol2inv(chol(K)) %>% as("dsCMatrix"))
        temp <- cholPermute(Qx)
        ybar <- t(Sm@S) %*% Qobs %*% (Sm@Z - CZ %*% X %*% alpha)
        x_mean <- cholsolve(Qx,ybar,perm=TRUE,cholQp = temp$Qpermchol, P = temp$P)
        if(predict_BAUs) {
            Cov <- Takahashi_Davis(Qx,cholQp = temp$Qpermchol,P = temp$P) # PARTIAL
        } else {
            Cov <- cholsolve(Qx,Diagonal(nrow(Qx)),perm=TRUE,
                             cholQp = temp$Qpermchol, P = temp$P) # FULL
        }
        BAUs[["mu"]] <- as.numeric(X %*% alpha) + as.numeric(S0 %*% x_mean)
        BAUs[["var"]] <- .batch_compute_var(S0,Cov,obs_fs = TRUE)
        #BAUs[["var"]] <- rowSums((S0 %*% Cov) * S0) + Sm@sigma2fshat*BAUs$fs
    }

    if(is.null(pred_polys)) {
        BAUs[["sd"]] <- sqrt(BAUs[["var"]])
        if(!is.null(pred_time)) {
            BAUs[,pred_time]
        } else {
            BAUs
        }
    } else {
        pred_polys[["mu"]] <- as.numeric(CP %*% BAUs[["mu"]])
        if(!obs_fs) CPM <- CP %*% PI else CPM <- CP %*% S0
        if(sigma2fs == 0)  pred_polys[["var"]] <- diag2(CPM %*% Cov, t(CPM)) ## All Cov available
        else pred_polys[["var"]] <- diag2(CPM, cholsolve(Q=Qx,y=t(CPM), ## Cov not available
                                                         perm = TRUE,cholQp = chol_Qx$Qpermchol,P = chol_Qx$P))
        pred_polys[["sd"]] <- sqrt(pred_polys[["var"]])
        pred_polys
    }

}

SRE.simulate <- function(S,obs_fs) {

    ## Still in development:
    print("Normalising basis function evaluations at BAUs ...")
    xx <- sqrt(rowSums(S0* S0))
    xx <- xx + 1*(xx == 0) ## Where there are no basis functions do not divide zero by zero..
    S0 <- S0/xx
    print("Done ...")

}

############ DEPRECATED #####################

.SRE.Mstep.deprecated <- function(Sm,alpha_OLS = FALSE) {

    mu_eta <- Sm@mu_eta
    S_eta <- Sm@S_eta
    alpha_init <- Sm@alphahat
    sigma2fs_init <- Sm@sigma2fshat

    K <- S_eta + tcrossprod(mu_eta)
    alpha <- alpha_init
    sigma2fs <- sigma2fs_init
    converged <- FALSE

    # Deprecated:
    #   Omega_diag1 <- diag2(Sm@S %*% as(S_eta,"dgeMatrix"),t(Sm@S)) +
    #                 diag2(Sm@S %*% mu_eta %*% t(mu_eta), t(Sm@S))

    R_eta <- chol(S_eta + tcrossprod(mu_eta))
    S_R_eta <- Sm@S %*% t(R_eta)
    Omega_diag1 <- rowSums(S_R_eta^2)

    if(all((a <- diag(Sm@Ve)) == a[1]) &
       all((b <- diag(Sm@Vfs)) == b[1]) &
       all(rowSums(Sm@Vfs) == a[1]))    {
        homoscedastic <- TRUE
    } else {
        homoscedastic <- FALSE
    }

    while(!converged) {
        J <- function(sigma2fs) {
            if(sigma2fs < 0) {
                return(Inf)
            } else {
                D <- sigma2fs*Sm@Vfs + Sm@Ve
                Dinv <- chol2inv(chol(D))
                DinvV <- Dinv %*% Sm@Vfs
                -(-0.5*tr(DinvV) +
                      0.5*tr(DinvV %*% Dinv %*% Omega_diag)
                )
            }
        }

        resid <- Sm@Z - Sm@X %*% alpha
        Omega_diag <- Omega_diag1 -
            2*diag2(Sm@S %*% mu_eta, t(resid)) +
            diag2(resid,t(resid))
        Omega_diag <- Diagonal(x=Omega_diag)

        # Repeat until finding values on opposite sides of zero if heteroscedastic
        if(!homoscedastic) {
            amp_factor <- 10; OK <- 0
            while(!OK) {
                amp_factor <- amp_factor * 10
                if(!(sign(J(sigma2fs/amp_factor)) == sign(J(sigma2fs*amp_factor)))) OK <- 1
                if(amp_factor > 1e9) {
                    #warning("sigma2fs is being estimated to zero.
                    #        This might because because of an incorrect binning procedure.")
                    OK <- 1
                }
            }

            if(amp_factor > 1e9) {
                sigma2fs_new <- 0
                converged <- TRUE
            }
            sigma2fs_new <- stats::uniroot(f = J,
                                           interval = c(sigma2fs/amp_factor,sigma2fs*amp_factor))$root
        } else {
            sigma2fs_new <- 1/b[1]*(sum(Omega_diag)/length(Sm@Z) - a[1])
        }
        D <- sigma2fs_new*Sm@Vfs + Sm@Ve
        Dinv <- chol2inv(chol(D))
        if(alpha_OLS) {
            converged <- TRUE
        } else {
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



.batch_compute_var.deprecated <- function(X,Cov) {
    # Don't consider more than 50e6 elements at a time
    batch_size <- min(round(50e6 / nrow(Cov)),nrow(Cov))
    batching=cut(1:nrow(X),breaks = seq(0,nrow(X)+batch_size,by=batch_size),labels=F)
    if(opts_FRK$get("parallel") > 1 & batch_size < nrow(X)) {
        clusterExport(opts_FRK$get("cl"),
                      c("batching","X","Cov"),envir=environment())
        var_list <- parLapply(opts_FRK$get("cl"),1:max(unique(batching)),
                              function(i) {
                                  idx = which(batching == i)
                                  as.numeric(rowSums((X[idx,] %*% Cov)*X[idx,]))})
        clusterEvalQ(opts_FRK$get("cl"), {gc()})

        temp <- do.call(c,var_list)
    } else {
        temp <- rep(0,nrow(X))
        for(i in 1:max(unique(batching))) {
            idx = which(batching==i)
            temp[idx] <- as.numeric(rowSums((X[idx,] %*% Cov)*X[idx,]))
        }
    }
    temp
}

#' @title Convert data frame to SpatialPolygons
#' @description Convert data frame to SpatialPolygons object.
#' @param df data frame containing polygon information, see details
#' @param keys vector of variable names used to group rows belonging to the same polygon
#' @param coords vector of variable names identifying the coordinate columns
#' @param proj the projection of the \code{SpatialPolygons} object. Needs to be of class \code{CRS}
#' @details Each row in the data frame \code{df} contains both coordinates and labels (or keys) that identify to which polygon the coordinates belong. This function groups the data frame according to \code{keys} and forms a \code{SpatialPolygons} object from the coordinates in each group. It is important that all rings are closed, that is, that the last row of each group is identical to the first row. Since \code{keys} can be of length greater than one, we identify each polygon with a new key by forming an MD5 hash made out of the respective \code{keys} variables that in themselves are unique (and therefore the hashed key is also unique). For lon-lat coordinates use \code{proj = CRS("+proj=longlat +ellps=sphere")}.
#' @export
#' @examples
#' library(sp)
#' df <- data.frame(id = c(rep(1,4),rep(2,4)),
#'                  x = c(0,1,0,0,2,3,2,2),
#'                  y=c(0,0,1,0,0,1,1,0))
#' pols <- df_to_SpatialPolygons(df,"id",c("x","y"),CRS())
#' \dontrun{plot(pols)}
df_to_SpatialPolygons2 <- function(df,keys,coords,proj) {

    ## Basic checks
    if(!is(df,"data.frame")) stop("df needs to be a data frame")
    if(!is(keys,"character")) stop("keys needs to be of class character")
    if(!is(coords,"character")) stop("coords needs to be of class character")
    if(!all(keys %in% names(df))) stop("All keys needs to be labels in data frame")
    if(!all(coords %in% names(df))) stop("All coordinate labels needs to be labels in data frame")
    if(!is(proj,"CRS")) stop("proj needs to be of class CRS")

    ## dfun takes a data frame with coordinates for 1 polygon, and makes one POLYGON object from it
    ## with a UID from the polygon key
    dfun <- function(d) {
        Polygons(list(Polygon(d[coords])),digest::digest(d[keys]))
    }

    if(0) { ## Do not enable, mostly overhead
        ## Deprecated to remove plyr:
        #doParallel::registerDoParallel(opts_FRK$get("parallel"))
        #df_poly <- plyr::dlply(df,keys,dfun,.parallel=TRUE)

        unique_keys <- unique(data.frame(df[keys]))[,1]

        clusterExport(opts_FRK$get("cl"),
                      c("df"),envir=environment())
        df_poly <- parLapply(opts_FRK$get("cl"),unique_keys,
                             function(key) {
                                 df[df[keys]==key,] %>%
                                     data.frame() %>%
                                     dfun})
        clusterEvalQ(opts_FRK$get("cl"), {gc()})

        # df_poly <- mclapply(unique_keys,
        #                     function(key) {
        #                         df[df[keys]==key,] %>%
        #                             data.frame() %>%
        #                             dfun},
        #                     mc.cores = opts_FRK$get("parallel"))

    } else {
        df_poly <- plyr::dlply(df,keys,dfun)
    }

    ## Rhipe version (currently disabled)

    # df_poly <- rhwrapper(Ntot = nrow(df),
    #                      N = 4000,
    #                      f_expr = .rhdlply,
    #                      df=df,
    #                      keys=keys,
    #                      coords=coords,
    #                      dfun=parse(text = deparse(dfun)))

    Sr <- SpatialPolygons(df_poly,1:length(df_poly),proj4string=proj)
}

.prec_from_neighb <- function (neighb, intrinsic = 1, precinc = 1)
{
    num_v <- length(neighb)
    num_neighb <- lapply(neighb, length)
    if (intrinsic == 1) {
        i_list <- vector("list", num_v)
        for (k in 1:num_v) {
            i_list[[k]] <- rep(k, num_neighb[[k]])
        }
        i <- unlist(i_list)
        j <- unlist(neighb)
        z <- rep(-1, length(j))
        i <- c(i, 1:num_v)
        j <- c(j, 1:num_v)
        zdiag <- unlist(num_neighb)
        z <- c(z, zdiag)
    }
    if (intrinsic == 2) {
        i1 <- 1:num_v
        j1 <- 1:num_v
        z1 <- rep(0, num_v)
        for (k in 1:num_v) {
            z1[k] <- num_neighb[[k]]^2 + num_neighb[[k]]
        }
        count <- 1
        i2 <- rep(0, num_v * 10)
        j2 <- rep(0, num_v * 10)
        z2 <- rep(0, num_v * 10)
        for (k in 1:num_v) {
            for (l in neighb[[k]]) {
                i2[count] <- k
                j2[count] <- l
                z2[count] <- -(num_neighb[[k]] + num_neighb[[l]] -
                                   sum(duplicated(c(neighb[[l]], neighb[[k]]))))
                count <- count + 1
            }
        }
        i2 <- i2[1:count - 1]
        j2 <- j2[1:count - 1]
        z2 <- z2[1:count - 1]
        count <- 1
        i3 <- rep(0, num_v * 15)
        j3 <- rep(0, num_v * 15)
        z3 <- rep(0, num_v * 15)
        neighb2 <- vector("list", num_v)
        for (k in 1:num_v) {
            for (l in neighb[[k]]) {
                neighb2[[k]] <- c(neighb2[[k]], setdiff(neighb[[l]],
                                                        c(neighb[[k]], k)))
            }
            for (l in unique(neighb2[[k]])) {
                i3[count] <- k
                j3[count] <- l
                z3[count] <- sum(neighb2[[k]] == l)
                count <- count + 1
            }
        }
        i3 <- i3[1:count - 1]
        j3 <- j3[1:count - 1]
        z3 <- z3[1:count - 1]
        i <- c(i1, i2, i3)
        j <- c(j1, j2, j3)
        z <- c(z1, z2, z3)
    }
    z <- precinc * z
    Q <- sparseMatrix(i, j, x = z)
    return(Q)
}

