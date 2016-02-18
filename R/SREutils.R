#' @title Construct SRE object, fit and predict
#' @description Main constructor of spatial random effects (SRE) object. Please see \code{\link{SRE-class}} for more details on the object's properties and methods.
#' @param f \code{R} formula relating the dependent variable (or transformations thereof) to covariates
#' @param data list of objects of class \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}
#' @param basis object of class \code{Basis}
#' @param BAUs object of class \code{SpatialPolygonsDataFrame}, the data frame which must contain covariate information as well as a field \code{fs} describing the fine-scale variation up to a constant of proportionality
#' @param est_error flag indicating whether the measurement-error variance should be estimated from variogram techniques. If this is set to 0, then \code{data} must contain a field \code{std}. Measurement-error estimation is not implemented for spatio-temporal datasets
#' @param average_in_BAU if \code{TRUE}, then multiple data points falling in the same BAU are averaged; the measurement error of the averaged data point is taken as the average of the individual measurement errors
#' @param SRE_model object returned from the constructor \code{SRE()}
#' @param n_EM maximum number of iterations for the EM algorithm
#' @param tol convergence tolerance for the EM algorithm
#' @param method parameter estimation method to employ. Currently only ``EM'' is supported
#' @param print_lik flag indicating whether likelihood should be printed or not on convergence of the estimation algorithm
#' @param use_centroid flag indicating whether the basis functions are averaged over the BAU, or whether the basis functions are evaluated at the BAUs centroid in order to construct the matrix \eqn{S}. The flag can safely be set when the basis functions are approximately constant over the BAUs in order to reduce computational time
#' @param obs_fs flag indicating whether the fine-scale variation sits in the observation model (systematic error) or in the process model (process fine-scale variation)
#' @param pred_polys object of class \code{SpatialPoylgons} indicating the regions over which prediction will be carried out. The BAUs are used if this option is not specified
#' @param pred_time vector of time indices at which we wish to predict. All time points are used if this option is not specified
#' @details \code{SRE()} is the main function in the package as it constructs a spatial random effects model from the user-defined formula, data object, basis functions and a set of Basic Areal Units (BAUs). The function first takes each object in the list \code{data} and maps it to the BAUs -- this entails binning the point-referenced data into BAUs (and averaging within the BAU) if \code{average_in_BAU = TRUE}, and finding which BAUs are influenced by the polygon datasets. Following this, the incidence matrix \code{Cmat} is constructed, which appears in the observation model \eqn{Z = CY + e}, where \eqn{C} is the incidence matrix.
#'
#' The SRE model is given by \eqn{Y = T\alpha + S\eta + \delta}, where \eqn{X} are the covariates at BAU level, \eqn{\alpha} are the regression coefficients, \eqn{S} are the basis functions evaluated at the BAU level, \eqn{\eta} are the basis function weights, and \eqn{\delta} is the fine scale variation (at the BAU level). The covariance matrix of \eqn{\delta} is diagonal and proportional to the field `fs' in the BAUs (typically set to one). The constant of proportionality is estimated in the EM algorithm. All required matrices (\eqn{S,T} etc.) are computed and returned as part of the object, please see \code{\link{SRE-class}} for more details.
#'
#'\code{SRE.fit()} takes an object of class \code{SRE} and estimates all unknown parameters, namely the covariance matrix \eqn{K}, the fine scale variance \eqn{\sigma^2_{\delta}} and the regression parameters \eqn{\alpha}. The only method currently implemented is the Expectation Maximisation (EM) algorithm, which the user configures through \code{n_EM} and \code{tol}. The latter parameter, \code{tol}, is used as in Katzfuss and Cressie to, that is, the log-likelihood (given in Equation (16) in that work) is evaluated at each iteration at the current parameter estimate, and convergence is assumed to have been reached when this quantity stops changing by more than \code{tol}.
#'
#'The actual computations for the E-step and M-step are relatively straightforward. The E-step contains an inverse of an \eqn{n \times n} matrix, where \code{n} is the number of basis functions which should not exceed 2000. The M-step first updates the matrix \eqn{K}, which only depends on the sufficient statistics of the basis weights \eqn{\eta}. Then, the regression parameter \eqn{\alpha} is updated and a simple optimisation routine (a line search) is used to update the fine-scale variance \eqn{\sigma^2_{\delta}}. If the fine-scale errors and measurement errors are homoscedastic a closed-form solution is available for the update of \eqn{\sigma^2_{fs}}. Irrespectively, since the udpates of \eqn{\alpha} and \eqn{\sigma^2_{\delta}} are dependent, these two updates are iterated until the change in \eqn{\sigma^2_{\delta}} is no more than 0.1\%.
#'
#'Once the parameters are fitted, the \code{SRE} object is passed onto the function \code{SRE.predict()} in order to carry out optimal predictions over the same BAUs used to construct the SRE model with \code{SRE()}. The first part of the prediction process is to construct the matrix \eqn{S} by averaging the basis functions over the prediction polygons/BAUs, using Monte Carlo integration with 1000 samples. This is a computationally-intensive process. On the other hand, one could set \code{use_centroid = TRUE} and treat the prediction over the entire polygon as that of a BAU at the centre of the polygon. This will yield valid results only if the BAUs are relatively small. Once the matrix \eqn{S} is found, a standard Gaussian inversion using the estimated parameters is used.
#'
#'\code{SRE.predict} returns the BAUs, which are of class \code{SpatialPolygonsDataFrame}, with two added attributes, \code{mu} and \code{var}. These can then be easily plotted using \code{spplot} or \code{ggplot2} (in conjunction with \code{\link{SpatialPolygonsDataFrame_to_df}}) as shown in the package vignettes.
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
#' grid_BAUs <- auto_BAUs(manifold=real_line(),data=sim_data,cellsize = c(0.01),type="grid")
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
#' S <- SRE.fit(S,n_EM = 5,tol = 1e-5,print_lik=TRUE)
#'
#' ### Predict over BAUs
#' grid_BAUs <- SRE.predict(S,use_centroid = TRUE)
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
SRE <- function(f,data,basis,BAUs,est_error=FALSE,average_in_BAU = TRUE) {

    .check_args(f=f,data=data,basis=basis,BAUs=BAUs,est_error=est_error)
    av_var <-all.vars(f)[1]
    ndata <- length(data)

    S <- Ve <- Vfs <- X <- Z <- Cmat <- list()
    for(i in 1:ndata) {
        if(est_error) data[[i]]$std <- 0 ## Just set it to something, this will be overwritten later on

        data_proc <- map_data_to_BAUs(data[[i]],
                                      BAUs,
                                      av_var = av_var,
                                      average_in_BAU = average_in_BAU)

        if(any(is.na(data_proc@data[av_var])))
            stop("NAs found when mapping data to BAUs. Are you sure all your data are covered by BAUs?")

        if(est_error) {
            if(is(data_proc,"ST"))
                stop("Estimation of error not yet implemented for spatio-temporal data")
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

        if(any(rowSums(Cmat[[i]])==0))
            stop("I have found difficulty in associating the data with the BAUs. If you have point-referenced data
                 then this could be because you have data outside BAUs. If you have polygon data, then
                 this could be because no BAUs centroids are within the polygons. For polygon data,
                 influence on a BAU is determined from whether the BAU centroid falls within the polygon
                 or not.")

        Cmat[[i]] <- Cmat[[i]] / rowSums(Cmat[[i]]) ## Average BAUs for polygon observations

        Vfs[[i]] <- tcrossprod(Cmat[[i]] %*% Diagonal(x=sqrt(BAUs$fs)))
        # if(length(Cmat[[i]]@x) == nrow(Cmat[[i]]))  { # if point observations
        #     Vfs[[i]] <- Diagonal(x = Cmat[[i]] %*% sqrt(BAUs$fs) %>% as.numeric())
        # } else {
        #     Vfs[[i]] <- Diagonal(x = rep(0,nrow(Cmat[[i]])))
        # }

        print("Evaluating basis functions at observation locations...")
        S[[i]] <- eval_basis(basis, s = data_proc)
        print("Done.")

        ## Still in development:
        # print("Normalising basis function evaluations at BAUs ...")
        # S0 <- eval_basis(basis,as.matrix(BAUs[coordnames(Sm@data[[1]])]@data))
        # xx <- sqrt(rowSums((S0) * S0))
        # S[[i]] <- S[[i]] / (Cmat[[i]] %*% xx %>% as.numeric())
        # print("Done ...")

        ## Note that S constructed in this way is similar to Cmat %*% S_BAUs where S_BAUs is the
        ## basis functions evaluated at the BAUs. Verify this by checking the following are similar
        ## S2 <- eval_basis(basis, s = BAUs)
        ## (Cmat[[i]] %*% S2) - S[[i]]
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
        Q_eta = Diagonal(x = rep(1,nbasis(basis))),
        alphahat = solve(t(X) %*% X) %*% t(X) %*% Z,
        Khat = Diagonal(n=nbasis(basis),x = var(Z[,1])),
        #Khat_inv = Diagonal(n=nbasis(basis),x = 1/var(Z[,1])),
        #B_run = Diagonal(n=nbasis(basis),x = 1/var(Z[,1])),
        #v_run = Matrix(0,nbasis(basis),nbasis(basis)),
        sigma2fshat = mean(diag(Ve)) / mean(diag(Vfs)))
}

#' @rdname SRE
#' @export
SRE.fit <- function(SRE_model,n_EM = 100L, tol = 1e-5, method="EM", print_lik=FALSE) {
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
        plot(1:i,lk[1:i],ylab="log likelihood",xlab="EM iteration (from #1)")
    }
    SRE_model
}


.SRE.Estep <- function(Sm) {

    alpha <- Sm@alphahat
    K <- Sm@Khat
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


    Kinv <- chol2inv(chol(K))

    # Deprecated: S_eta <- chol2inv(chol(crossprod(sqrt(Dinv) %*% Sm@S) + Kinv))
    S_eta <- chol2inv(chol(crossprod(t(cholDinv) %*% Sm@S) + Kinv))
    mu_eta <- S_eta %*% (t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha))


    # Sm@B_run <- Sm@B_run +  crossprod(t(cholDinv) %*% Sm@S)
    # Q_eta2 <- crossprod(t(cholDinv) %*% Sm@S) + Kinv
    # Q_eta <- Sm@B_run - Sm@v_run
    # mu_eta2 <- solve(Q_eta,(t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha)))
    # muQ <- Q_eta %*% mu_eta
    # c <- (1 + t(mu_eta) %*% muQ) %>% as.numeric()
    # v <- tcrossprod(muQ,muQ)
    # Sm@v_run <- Sm@v_run + 1/c * v



    Sm@mu_eta <- mu_eta
    Sm@S_eta <- S_eta
    # Sm@Q_eta <- Q_eta

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

    if(all((a <- diag(Sm@Ve)) == a[1]) &
       all((b <- diag(Sm@Vfs)) == b[1]) &
       all(rowSums(Sm@Vfs) == b[1]))    {
        homoscedastic <- TRUE
    } else {
        homoscedastic <- FALSE
    }

    if( all(rowSums(Sm@Ve) == diag(Sm@Ve)) &
        all(rowSums(Sm@Vfs) == diag(Sm@Vfs)))    {
        diagonal_mats <- TRUE
    } else {
        diagonal_mats <- FALSE
    }



    if(!diagonal_mats) {

        J <- function(sigma2fs) {
            if(sigma2fs < 0) {
                return(Inf)
            } else {
                D <- sigma2fs*Sm@Vfs + Sm@Ve
                Dinv <- chol2inv(chol(D))
                DinvV <- Dinv %*% Sm@Vfs
                DinvVDinv <- Dinv %*% Sm@Vfs %*% Dinv

                alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*% t(Sm@X) %*% Dinv %*% (Sm@Z - Sm@S %*% mu_eta)
                resid <- Sm@Z - Sm@X %*% alpha

                Dinvr <- DinvVDinv %*% resid
                DinvS <- DinvVDinv %*% Sm@S

                tr1 <- tr(DinvV)
                tr2 <- sum(diag2(DinvS %*% (S_eta +  tcrossprod(mu_eta)),t(Sm@S))  -
                               2*diag2(DinvS %*% mu_eta,t(resid)) +
                               diag2(Dinvr,t(resid)))

                -(-0.5*tr1 +0.5*tr2)
            }
        }
    } else {

        R_eta <- chol(S_eta + tcrossprod(mu_eta))
        S_R_eta <- Sm@S %*% t(R_eta)
        Omega_diag1 <- rowSums(S_R_eta^2)

        J <- function(sigma2fs) {
            if(sigma2fs < 0) {
                return(Inf)
            } else {
                D <- sigma2fs*Sm@Vfs + Sm@Ve
                if(is(D,"dtCMatrix")) {
                    Dinv <- solve(D)
                } else {
                    Dinv <- chol2inv(chol(D))
                }
                DinvV <- Dinv %*% Sm@Vfs

                alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*% t(Sm@X) %*% Dinv %*% (Sm@Z - Sm@S %*% mu_eta)
                resid <- Sm@Z - Sm@X %*% alpha
                Omega_diag <- Omega_diag1 -
                    2*diag2(Sm@S %*% mu_eta, t(resid)) +
                    diag2(resid,t(resid))
                Omega_diag <- Diagonal(x=Omega_diag)

                -(-0.5*tr(DinvV) +
                  0.5*tr(DinvV %*% Dinv %*% Omega_diag)
                )
            }
        }
    }



    # Repeat until finding values on opposite sides of zero if heteroscedastic
    if(!homoscedastic) {
        amp_factor <- 10; OK <- 0
        while(!OK) {
            amp_factor <- amp_factor * 10
            if(!(sign(J(sigma2fs/amp_factor)) == sign(J(sigma2fs*amp_factor)))) OK <- 1
            if(amp_factor > 1e12) {
                warning("sigma2fs is being estimated to zero.
                            This might because because of an incorrect binning procedure.")
                OK <- 1
            }
        }

        if(amp_factor > 1e12) {
            sigma2fs_new <- 0
        }
        sigma2fs_new <- uniroot(f = J,
                                interval = c(sigma2fs/amp_factor,sigma2fs*amp_factor))$root
        D <- sigma2fs_new*Sm@Vfs + Sm@Ve
        if(is(D,"dtCMatrix")) {
            Dinv <- solve(D)
        } else {
            Dinv <- chol2inv(chol(D))
        }
        alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*% t(Sm@X) %*% Dinv %*% (Sm@Z - Sm@S %*% mu_eta)
    } else {

        alpha <- solve(t(Sm@X) %*% Sm@X) %*% t(Sm@X) %*% (Sm@Z - Sm@S %*% mu_eta)
        resid <- Sm@Z - Sm@X %*% alpha
        Omega_diag <- Omega_diag1 -
            2*diag2(Sm@S %*% mu_eta, t(resid)) +
            diag2(resid,t(resid))
        Omega_diag <- Diagonal(x=Omega_diag)
        sigma2fs_new <- 1/b[1]*(sum(Omega_diag)/length(Sm@Z) - a[1])
        if(sigma2fs_new < 0) {
            warning("sigma2fs is being estimated to zero.
                            This might because because of an incorrect binning procedure or
                    because too much measurement error is being assumed.")
            sigma2fs_new = 0
        }
    }

    Sm@Khat <- K
    Sm@alphahat <- alpha
    Sm@sigma2fshat <- sigma2fs_new

    Sm
}

.loglik <- function(Sm) {

    D <- Sm@sigma2fshat*Sm@Vfs + Sm@Ve
    S <- Sm@S
    K <- Sm@Khat
    chol_K <- chol(K)
    Kinv <- chol2inv(chol_K)
    resid <- Sm@Z - Sm@X %*% Sm@alphahat
    N <- length(Sm@Z)

    if(is(D,"dtCMatrix")) {
        cholD <- sqrt(D)
        cholDinvT <- solve(cholD)
    } else {
        cholD <- chol(D)
        cholDinvT <- t(solve(cholD))
    }


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


#' @rdname SRE
#' @export
SRE.predict <- function(SRE_model,use_centroid=TRUE,obs_fs=TRUE,pred_polys = NULL,pred_time = NULL) {


     if(is(pred_polys,"Spatial") & !(is(pred_polys,"SpatialPolygons")))
         stop("Predictions need to be over BAUs or spatial polygons")

     if(is(pred_polys,"ST") & !(is(pred_polys,"STFDF")))
         if(!(is(pred_polys@sp,"SpatialPolygons")))
             stop("Predictions need to be over BAUs or STFDFs with spatial polygons")

    pred_locs <- .SRE.predict(Sm=SRE_model,
                              use_centroid=use_centroid,
                              obs_fs=obs_fs,
                              pred_polys = pred_polys,
                              pred_time = pred_time)

    ## Rhipe VERSION (currently disabled)
    # pred_locs <- rhwrapper(Ntot = length(Sm@BAUs),
    #                        N = 4000,
    #                        f_expr = .rhSRE.predict,
    #                        Sm = SRE_model,
    #                        pred_locs = pred_locs,
    #                        use_centroid = use_centroid)

    pred_locs
}

.SRE.predict <- function(Sm,use_centroid,obs_fs = TRUE,pred_polys = NULL,pred_time = NULL) {

    if(!is.null(pred_polys))
        if(is(Sm@BAUs,"ST"))
            stop("Prediciton is currently only possible at BAU level with
                  spatio-temporal models. Please do not use pred_polys")

    if(is.null(pred_time) & is(Sm@BAUs,"ST"))
        pred_time <- 1:length(Sm@BAUs@time)


    predict_BAUs <- TRUE
    BAUs <- Sm@BAUs

    if(is.null(pred_polys)) {
        CP <- Diagonal(length(BAUs))
    } else {
        C_idx <- BuildC(pred_polys,BAUs)
        CP <- sparseMatrix(i=C_idx$i_idx,
                           j=C_idx$j_idx,
                           x=1,
                           dims=c(length(pred_polys),
                                      length(BAUs)))
        CP <- CP / rowSums(CP) ## Average over polygon
        if(!all(table(C_idx$i_idx) == 1))
            predict_BAUs <- FALSE   ## Need to compute full covariance matrix
    }

    CZ <- Sm@Cmat

    ## Now, we only need those BAUs that are influenced by observations and prediction locations
    if(!predict_BAUs) {
        needed_BAUs <- union(as(CP,"dgTMatrix")@j+1, as(CZ,"dgTMatrix")@j+1)
        BAUs <- BAUs[needed_BAUs,]
        CP <- CP[,needed_BAUs]
        CZ <- CZ[,needed_BAUs]
    }

    # if(is(BAUs,"ST")){
    #     needed_BAUs <- BAUs[,pred_time]$n
    #     BAUs <- BAUs[,pred_time]
    #     CP <- CP[,needed_BAUs]
    #     CZ <- CZ[,needed_BAUs]
    # }

    depname <- all.vars(Sm@f)[1]
    BAUs[[depname]] <- 0.1
    L <- .gstat.formula(Sm@f,data=BAUs)
    X = as(L$X,"Matrix")
    if(is(BAUs,"Spatial")) {
        if(use_centroid) {
            S0 <- eval_basis(Sm@basis,as.matrix(BAUs[coordnames(Sm@data[[1]])]@data))
        } else {
            S0 <- eval_basis(Sm@basis,BAUs)
        }
    } else if(is(BAUs,"STFDF")) {
        if(use_centroid) {
            S0 <- eval_basis(Sm@basis,as.matrix(cbind(coordinates(BAUs),BAUs@data$t)))
        } else {
            stop("Can only use centroid when predicting with spatio-temporal data")
        }
    }

    ## Still in development:
    # print("Normalising basis function evaluations at BAUs ...")
    # xx <- sqrt(rowSums((S0) * S0))
    # S0 <- S0/xx
    # print("Done ...")

    BAUs[[depname]] <- NULL

    alpha <- Sm@alphahat
    K <- Sm@Khat
    sigma2fs <- Sm@sigma2fshat
    D <- sigma2fs*Sm@Vfs + Sm@Ve
    if(is(D,"dtCMatrix")) {
        Dchol <- sqrt(D)
        Dinv <- solve(D)
    } else {
        Dchol <- chol(D)
        Dinv <- chol2inv(Dchol)
    }
    mu_eta <- Sm@mu_eta
    S_eta <- Sm@S_eta

    sig2_Vfs_pred <- Diagonal(x=sigma2fs*BAUs$fs)

    if(is(BAUs,"Spatial")) {
        idx <- match(row.names(BAUs),row.names(Sm@BAUs))
    } else if (is(BAUs,"STFDF")){
        idx <- match(BAUs@data$n,Sm@BAUs@data$n)
    }



    if(!obs_fs) {
        if(sigma2fs >0) {
            LAMBDA <- as(bdiag(Sm@Khat,sig2_Vfs_pred),"symmetricMatrix")
            LAMBDAinv <- chol2inv(chol(LAMBDA))  #block diagonal so straightforward
            PI <- cBind(S0, .symDiagonal(n=length(BAUs)))
            Qx <- t(PI) %*% t(CZ) %*% solve(Sm@Ve) %*% CZ %*% PI + LAMBDAinv
            temp <- cholPermute(Qx)
            ybar <- t(PI) %*%t(CZ) %*% solve(Sm@Ve) %*% (Sm@Z - CZ %*% X %*% alpha)
            x_mean <- cholsolve(Qx,ybar,perm=TRUE,cholQp = temp$Qpermchol, P = temp$P)
            if(predict_BAUs) {
                Cov <- Takahashi_Davis(Qx,cholQp = temp$Qpermchol,P = temp$P) # PARTIAL
            } else {
                Cov <- cholsolve(Qx,Diagonal(nrow(Qx)),perm=TRUE,
                                 cholQp = temp$Qpermchol, P = temp$P) # FULL
            }
            BAUs[["mu"]] <- as.numeric(X %*% alpha + PI %*% x_mean)
        } else {
            LAMBDA <- as(Sm@Khat,"symmetricMatrix")
            LAMBDAinv <- chol2inv(chol(LAMBDA))
            PI <- S0
            Qx <- crossprod(solve(sqrt(Sm@Ve)) %*% CZ %*% PI) + LAMBDAinv
            ybar <- t(PI) %*%t(CZ) %*% solve(Sm@Ve) %*% (Sm@Z - CZ %*% X %*% alpha)
            Cov <- as(chol2inv(chol(Qx)),"dgeMatrix")  ## Do all covariance matrix
            ## We can do all the covariance matrix since the dimension is equal to those of eta
            x_mean <- Cov %*% ybar
            BAUs[["mu"]] <- as.numeric(X %*% alpha + PI %*% x_mean)
        }

        ## variance too hard to compute all at once -- do it in blocks of 1000
        batching=cut(1:nrow(PI),breaks = seq(0,nrow(PI)+1000,by=1000),labels=F)

        if(opts_FRK$get("parallel") > 1) {
            cl <- makeCluster(opts_FRK$get("parallel"))
            var_list <- mclapply(1:max(unique(batching)),
                                function(i) {
                                    idx = which(batching == i)
                                    as.numeric(rowSums((PI[idx,] %*% Cov)*PI[idx,]))},
                                mc.cores = opts_FRK$get("parallel"))
            temp <- do.call(c,var_list)
            stopCluster(cl)
        } else {
            temp <- rep(0,length(BAUs))
            for(i in 1:max(unique(batching))) {
                idx = which(batching==i)
                temp[idx] <- as.numeric(rowSums((PI[idx,] %*% Cov)*PI[idx,]))
            }
        }
        BAUs[["var"]] <- temp
    }

    if(obs_fs) {

        Qx <- (crossprod(t(solve(Dchol)) %*% (Sm@S %>% as("dgCMatrix"))) + chol2inv(chol(K)) %>% as("dsCMatrix"))
        temp <- cholPermute(Qx)
        ybar <- t(Sm@S) %*% Dinv %*% (Sm@Z - CZ %*% X %*% alpha)
        x_mean <- cholsolve(Qx,ybar,perm=TRUE,cholQp = temp$Qpermchol, P = temp$P)
        if(predict_BAUs) {
            Cov <- Takahashi_Davis(Qx,cholQp = temp$Qpermchol,P = temp$P) # PARTIAL
        } else {
            Cov <- cholsolve(Qx,Diagonal(nrow(Qx)),perm=TRUE,
                             cholQp = temp$Qpermchol, P = temp$P) # FULL
        }
        BAUs[["mu"]] <- as.numeric(X %*% alpha + S0 %*% x_mean)
        BAUs[["var"]] <- rowSums((S0 %*% Cov) * S0)
    }

    if(predict_BAUs) {
        if(!is.null(pred_time)) {
            BAUs[,pred_time]
        } else {
            BAUs
        }
    } else {
        pred_polys[["mu"]] <- as.numeric(CP %*% BAUs[["mu"]])
        if(!obs_fs) {
            pred_polys[["var"]] <- diag2(CP %*% PI %*% Cov, t(PI) %*% t(CP))
        } else {
            pred_polys[["var"]] <- diag2(CP %*% S0 %*% Cov, t(S0) %*% t(CP))
        }
        pred_polys
    }

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
              cat(paste0("Fine-scale variance proportionality constant [extract using object@sigma2fshat]: ",object@sigma2fshat))
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




.check_args <- function(f,data,basis,BAUs,est_error) {
    if(!is(f,"formula")) stop("f needs to be a formula.")
    if(is(BAUs,"Spatial"))
        if(!all(all.vars(f)[-1] %in% names(BAUs@data)))
            stop("All covariates need to be in the SpatialPolygons BAU object")
    if(is(BAUs,"ST"))
        if(!all(all.vars(f)[-1] %in% c(names(BAUs@data),coordnames(BAUs))))
            stop("All covariates need to be in the SpatialPolygons BAU object")
    if(!is(data,"list"))
        stop("Please supply a list of Spatial objects.")
    if(!all(sapply(data,function(x) is(x,"Spatial") | is(x,"ST"))))
        stop("All data list elements need to be of class Spatial or ST")
    if(!all(sapply(data,function(x) all.vars(f)[1] %in% names(x@data))))
        stop("All data list elements to have values for the dependent variable")
    if(!all(sapply(data,function(x) identical(proj4string(x), proj4string(BAUs)))))
        stop("Please ensure all data items and BAUs have the same coordinate reference system")
    if(!(is(basis,"Basis") | is(basis,"TensorP_Basis")))
        stop("basis needs to be of class Basis  or TensorP_Basis (package FRK)")
    if(!("fs" %in% names(BAUs@data))) {
        warning("BAUs should contain a field 'fs' containing a basis
                function for fine-scale variation. Setting basis function equal to one everywhere.")
        BAUs$fs <- 1
    }
    if(!(all(BAUs$fs >= 0)))
        stop("fine-scale variation basis function needs to be nonnegative everywhere")
    if(!(is(BAUs,"SpatialPolygonsDataFrame") | is(BAUs,"STFDF")))
        stop("BAUs should be a SpatialPolygonsDataFrame or a STFDF object")
    if(is(BAUs,"STFDF")) if(!is(BAUs@sp,"SpatialPolygonsDataFrame"))
        stop("The spatial component of the BAUs should be a SpatialPolygonsDataFrame")
    if((is(manifold(basis),"sphere")) & !all((coordnames(BAUs) == c("lon","lat"))))
        stop("Since a sphere is being used, please ensure that
             all coordinates (including those of BAUs) are in (lon,lat)")
    if(!est_error & !all(sapply(data,function(x) "std" %in% names(x@data))))
        stop("If observational error is not going to be estimated,
             please supply a field 'std' in the data objects")
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
                if(amp_factor > 1e12) {
                    warning("sigma2fs is being estimated to zero.
                            This might because because of an incorrect binning procedure.")
                    OK <- 1
                }
            }

            if(amp_factor > 1e12) {
                sigma2fs_new <- 0
                converged <- TRUE
            }
            sigma2fs_new <- uniroot(f = J,
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
