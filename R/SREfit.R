#' @rdname SRE
#' @export
SRE.fit <- function(object, n_EM = 100L, tol = 0.01, method = c("EM", "TMB"),
                    lambda = 0, print_lik = FALSE, optimiser = nlminb, 
                    known_sigma2fs = NULL, taper = NULL, 
                    simple_kriging_fixed = TRUE, ...) {
  
  ## Deprecation coercion
  tmp <- list(...)
  if(!is.null(tmp$SRE_model)) {
    object <- tmp$SRE_model
    warning("The argument 'SRE_model' is deprecated: Please use 'object'")
  }

    method <- match.arg(method)
    optimiser <- match.fun(optimiser)
    
    if (!is.null(known_sigma2fs)) object@sigma2fshat <- known_sigma2fs
    
    if (method == "TMB" & object@K_type == "block-exponential") {
      answer <- user_decision("You have selected method = 'TMB' and K_type = 'block-exponential'. While this combination is allowed, it is significantly more computationally demanding than K_type = 'precision'. Please enter Y if you would like to continue with the block-exponential formulation, or N if you would like to change to the sparse precision matrix formulation.\n")
      if (answer == "N") {
        object@K_type <- "precision"
        cat("Setting K-type = 'precision'.\n")
      }
    }


    ## Check the arguments are OK
    ## (Note that we cannot pass the SRE model, because it complicates things 
    ## for the FRK() wrapper)
    .check_args2(n_EM = n_EM, tol = tol, lambda = lambda,
                 method = method, print_lik = print_lik, 
                 fs_by_spatial_BAU = object@fs_by_spatial_BAU, 
                 response = object@response, K_type = object@K_type, link = object@link, 
                 known_sigma2fs = known_sigma2fs,
                 BAUs = object@BAUs,
                 optimiser = optimiser, taper = taper, 
                 simple_kriging_fixed = simple_kriging_fixed, 
                 ...) # control parameters to optimiser() 
    
    object@simple_kriging_fixed <- simple_kriging_fixed
    

    ## Call internal fitting function with checked arguments
    object <- .SRE.fit(object = object, n_EM = n_EM, tol = tol, 
                          method = method, lambda = lambda, print_lik = print_lik, 
                          optimiser = optimiser, known_sigma2fs = known_sigma2fs, taper = taper, ...)
    return(object)
}




# ---- NOT EXPORTED ----

.SRE.fit <- function(object, n_EM, tol, method, lambda, 
                     print_lik, optimiser, known_sigma2fs, taper, ...) {
  
  if(method == "EM") {
    object <- .EM_fit(object = object, n_EM = n_EM, lambda = lambda, 
                         tol = tol, print_lik = print_lik, known_sigma2fs = known_sigma2fs, taper = taper)
  } else if (method == "TMB") {
    object <- .TMB_fit(object, optimiser = optimiser, known_sigma2fs = known_sigma2fs, taper = taper, ...)
  } else {
    stop("No other estimation method implemented yet. Please use method = 'EM' or method = 'TMB'.")
  }
  
  object@method <- method
  return(object) # return fitted SRE model
}

# ---- EM fitting functions ----

.EM_fit <- function(object, n_EM, lambda, tol, print_lik, known_sigma2fs, taper) {
  
  info_fit <- list()      
  
  info_fit$time <- system.time({
  
  n <- nbasis(object)  # number of basis functions
  X <- object@X        # covariates
  
  info_fit$method <- "EM"  # updated info_fit
  llk <- rep(0,n_EM)       # log-likelihood
  
  ## If user wishes to show progress show progress bar
  if(opts_FRK$get("progress"))
    pb <- utils::txtProgressBar(min = 0, max = n_EM, style = 3)
  
  ## If the user has requested covariance tapering, construct the taper matrix
  if (!is.null(taper)) {
    beta   <- .tapering_params(D_matrices = object@D_basis, taper = taper)
    T_beta <- .T_beta_taper_matrix(D_matrices = object@D_basis, beta = beta)
    info_fit$taper <- "NULL"
  } else {
    T_beta <- 1 
  }
  
  ## For each EM iteration step
  for(i in 1:n_EM) {
    llk[i] <- loglik(object)                          # compute the log-lik
    object <- .SRE.Estep(object)                   # compute E-step
    object <- .SRE.Mstep(object, lambda, known_sigma2fs, T_beta) # compute M-step
    if(opts_FRK$get("progress"))
      utils::setTxtProgressBar(pb, i)                    # update progress bar
    if(i>1)                                              # If we're not on first iteration
      if(abs(llk[i] - llk[i-1]) < tol) {                 # Compute change in log-lik
        cat("Minimum tolerance reached\n")               # and stop if less than tol
        break
      }
  }
  
  if(opts_FRK$get("progress")) close(pb)           # close progress bar
  info_fit$num_iterations <- i                     # update fit info
  
  
  ## If zero fine-scale variation detected just make sure user knows.
  ## This can be symptomatic of poor fitting
  if(object@sigma2fshat == 0) {
    info_fit$sigma2fshat_equal_0 <- 1
    if(opts_FRK$get("verbose") > 0)
      message("sigma2fs is being estimated to zero.
                        This might because of an incorrect binning
                        procedure or because too much measurement error
                        is being assumed (or because the latent
                        field is indeed that smooth, but unlikely).")
  } else {
    info_fit$sigma2fshat_equal_0 <- 0
  }
  
  ## If we have reached max. iterations, tell the user
  if(i == n_EM) {
    cat("Maximum EM iterations reached\n")
    info_fit$converged <- 0    # update info_fit
  } else {
    info_fit$converged <- 1   # update info_fit
  }
  
  ## Plot log-lik vs EM iteration plot
  info_fit$plot_lik <- list(x = 1:i, llk = llk[1:i],
                            ylab = "log likelihood",
                            xlab = "EM iteration")
  
  ## If user wants to see the log-lik vs EM iteration plot, plot it
  if(print_lik & !is.na(tol)) {
    plot(1:i, llk[1:i],
         ylab = "log likelihood",
         xlab = "EM iteration")
    
  }
  }) # system.time brackets

  object@info_fit <- info_fit
  
  return(object)
}



## E-Step
.SRE.Estep <- function(Sm) {
    # This is structured this way so that extra models for fs-variation
    # can be implemented later
    if(Sm@fs_model == "ind")
        Sm <- .SRE.Estep.ind(Sm)
    else stop("E-step only for independent fs-variation model currently implemented")
}

## Compute log-likelihood for independent fs-variation model
.loglik.ind <- function(Sm) {

    S <- Sm@S                               # basis-function matrix
    K <- Sm@Khat                            # random-effects cov. matrix
    chol_K <- chol(K)                       # its Cholesky
    Kinv <- chol2inv(chol_K)                # random-effects prec. matrix
    resid <- Sm@Z - Sm@X %*% Sm@alphahat    # residuals at fitted estimates
    N <- length(Sm@Z)                       # number of data points

    D <- Sm@sigma2fshat*Sm@Vfs + Sm@Ve      # total variance of data
    if(isDiagonal(D)) {                     # if this is diagonal
        D <- Diagonal(x = D@x)              # cast to diagonal matrix
        cholD <- sqrt(D)                    # just compute sqrt
        cholDinvT <- solve(cholD)           # and the inverse is just the reciprocal
    } else {
        cholD <- chol(D)                    # otherwise do the Cholesky
        cholDinvT <- t(solve(cholD))        # find the transposed (lower) inverse of the factor
    }

    ## Compute log-determinant. This is given by a formula in Section 2.2
    S_Dinv_S <-  crossprod(cholDinvT %*% S)
    R <- chol(Kinv + S_Dinv_S)
    log_det_SigmaZ <- logdet(R) +
        determinant(K,logarithm = TRUE)$modulus +
        logdet(cholD)  # this computes the log-determinant of a matrix from its Cholesky factor

    ## Alternatively: (slower but more direct)
    # Dinv <- chol2inv(chol(D))
    # SigmaZ_inv <- Dinv - Dinv %*% S %*% solve(Kinv + S_Dinv_S) %*% t(S) %*% Dinv
    # SigmaZ_inv2 <- Dinv - tcrossprod(Dinv %*% S %*% solve(R))

    ## Compute efficiently rDinv <- t(resid) %*% Dinv
    rDinv <- crossprod(cholDinvT %*% resid,cholDinvT)

    ## Compute the quadratic portion of the log-lik
    ## This is the same as quad_bit <- rDinv %*% resid - tcrossprod(rDinv %*% S %*% solve(R))
    ## but more efficient
    quad_bit <- crossprod(cholDinvT %*% resid) - tcrossprod(rDinv %*% S %*% solve(R))

    ## Now just add the bits together
    llik <- -0.5 * N * log(2*pi) -
        0.5 * log_det_SigmaZ -
        0.5 * quad_bit

    return(as.numeric(llik)) 

}

## M-step
.SRE.Mstep <- function(Sm, lambda = 0, known_sigma2fs, T_beta) {
    # This is structured this way so that extra models for fs-variation
    # can be implemented later
    if(Sm@fs_model == "ind")
        Sm <- .SRE.Mstep.ind(Sm, lambda = lambda, known_sigma2fs = known_sigma2fs, T_beta = T_beta)
    else stop("M-step only for independent fs-variation model currently implemented")
}

## E-step for independent fs-variation model
.SRE.Estep.ind <- function(Sm, known_sigma2fs) {
  
    alpha <- Sm@alphahat           # current regression coefficients estimates
    K <- Sm@Khat                   # current random effects covariance matrix estimate
    Kinv <- Sm@Khat_inv            # current random effects precision matrix estimate
    sigma2fs <- Sm@sigma2fshat     # current fs-variation factor estimate

    D <- sigma2fs*Sm@Vfs + Sm@Ve   # total variance-covariance of Z
    if(isDiagonal(D)) {            # if this is diagonal
        D <- Diagonal(x=D@x)       # cast to diagonal
        cholD <- sqrt(D)           # then the Cholesky is the sqrt
        cholDinv <- solve(cholD)   # the inverse Cholesky is the inverse-sqrt
        Dinv <- solve(D)           # the inverse is just reciprocal of diagonal elements
    } else {
        cholD <- Matrix::chol(D)   # if not diagonal then do Cholesky
        cholDinv <- solve(cholD)   # invery Cholesky factor
        Dinv <- chol2inv(cholD)    # find the inverse from the Cholesky factor
    }

    ## The below are simple Gaussian updating equations in Section 2.2 of the vignette
    Q_eta <- (crossprod(t(cholDinv) %*% Sm@S) + Kinv)
    S_eta <- chol2inv(chol(Q_eta))  # we can invert since we are low rank in FRK
    mu_eta <- (S_eta) %*%(t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha))

    Sm@mu_eta <- mu_eta  # update conditional mean
    Sm@S_eta <- S_eta    # update conditional covariance
    Sm@Q_eta <- Q_eta    # update conditional precision
    Sm                   # return SRE object
}

## M-step for the indepdent fine-scale variation model
.SRE.Mstep.ind <- function(Sm, lambda = 0, known_sigma2fs, T_beta) {

    mu_eta <- Sm@mu_eta              # current cond. mean of random effects
    S_eta <- Sm@S_eta                # current cond. cov. matrix of random effects
    alpha <- Sm@alphahat             # regression coefficients
    sigma2fs <- Sm@sigma2fshat       # fine-scale variance

    K <- .update_K(Sm,method=Sm@K_type,  # update the prior covariance matrix K
                   lambda = lambda)
    K <- K * T_beta # apply covariance taper, or just multiply by 1 if taper is not requested
    Khat_inv <- chol2inv(chol(K))        # compute the precision

    ## If the measurement and fs. variational covariance matricies
    ## are proportional to the identity then we have the
    ## special case of homoscedasticity
    homoscedastic <- all((a <- diag(Sm@Ve)) == a[1]) & 
      all((b <- diag(Sm@Vfs)) == b[1]) &
      isDiagonal(Sm@Vfs)   

    ## If the measurement and fs. variational covariance matricies
    ## are diagonal then we have another special case
    diagonal_mats <- isDiagonal(Sm@Ve) & isDiagonal(Sm@Vfs)  
    
    ## If the user has supplied a known value for the fine-scale variance, we 
    ## don't need to estimate sigma2fs.
    if(!is.null(known_sigma2fs)) {
        est_sigma2fs <- FALSE
        sigma2fs_new <- known_sigma2fs
    } else {
       est_sigma2fs <- TRUE # and sigma2fs_new will be estimated below 
    }
      
    ## If we have some fine-scale variation terms
    if(!all(diag(Sm@Vfs) == 0))
        ## And we're not in the diagonal case (this is the most comp. intensive)
        if(!diagonal_mats) {
            ## We first need to create a function whose root is sigma2fshat
            ## See 2.2 of vignette for equation details
            J <- function(sigma2fs) {
                if(sigma2fs < 0) {
                    return(Inf)                              # cannot be less than 0
                } else {
                    D <- sigma2fs*Sm@Vfs + Sm@Ve             # total data variance-covariance
                    Dinv <- chol2inv(chol(D))                # it's inverse (this will be blocked so still sparse)
                    DinvV <- Dinv %*% Sm@Vfs                 # summary matrix
                    DinvVDinv <- Dinv %*% Sm@Vfs %*% Dinv    # summary matrix

                    alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*%
                        t(Sm@X) %*% Dinv %*%                 # regression coefficients GLS estimates
                        (Sm@Z - Sm@S %*% mu_eta)
                    resid <- Sm@Z - Sm@X %*% alpha           # fitted residuals

                    Dinvr <- DinvVDinv %*% resid             # summary vector
                    DinvS <- DinvVDinv %*% Sm@S              # summary vector

                    tr1 <- tr(DinvV)                         # compute trace of first term

                    ## Compute trace of second term
                    tr2 <- sum(diag2(DinvS %*% (S_eta +  tcrossprod(mu_eta)),t(Sm@S))  -
                                   2*diag2(DinvS %*% mu_eta,t(resid)) +
                                   diag2(Dinvr,t(resid)))

                    ## return value of function
                    -(-0.5*tr1 +0.5*tr2)
                }
            }
        } else {

            ## If we have diagonal matrices then some simplifications are possible
            R_eta <- chol(S_eta + tcrossprod(mu_eta)) # Cholesky factor of E(eta eta^T)
            S_R_eta <- Sm@S %*% t(R_eta)              # summary matrix
            Omega_diag1 <- rowSums(S_R_eta^2)         # first part of diag(Omega) as in vignette
            J <- function(sigma2fs) {
                if(sigma2fs < 0) {
                    return(Inf)                       # sigma2fs >= 0
                } else {
                    D <- sigma2fs*Sm@Vfs + Sm@Ve    # total data variance. This must be diagonal
                    Dinv <- solve(D)                # just take reciprocal since D is defo. diagonal here
                    DinvV <- Dinv %*% Sm@Vfs        # summary matrix (diagonal)

                    alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*%
                        t(Sm@X) %*% Dinv %*%          # regression coefficients GLS estimates
                        (Sm@Z - Sm@S %*% mu_eta)
                    resid <- Sm@Z - Sm@X %*% alpha    # fitted residuals
                    Omega_diag <- Omega_diag1 -       # other parts of diag(Omega)
                        2*diag2(Sm@S %*% mu_eta, t(resid)) +
                        diag2(resid,t(resid))
                    Omega_diag <- Diagonal(x=Omega_diag)  # compute a diagonal Omega matrix

                    ## Since DinvV and Dinv are diagonal and we only want the trace,
                    ## we only need the diagonal elements of Omega in the following
                    return(-(-0.5*tr(DinvV) + 0.5*tr(DinvV %*% Dinv %*% Omega_diag)))
                }
            }
        }

    ## We need to find the root in J. For this we need to start uniroot with
    ## values on either side of sigma2fshat. The below implements
    ## a simple search algorithm for finding a good starting values.

    ## If we have some fine-scale variation terms
    if(!all(diag(Sm@Vfs) == 0)) {
        ## And we're not in the special homoscedastic case
        if(!homoscedastic) {
            if(est_sigma2fs) {
                amp_factor <- 10; OK <- 0  # initialise
                while(!OK) {
                    amp_factor <- amp_factor * 10 # widen the interval

                    ## If the signs are different, then we're OK, otherwise not
                    if(!(sign(J(sigma2fs/amp_factor)) == sign(J(sigma2fs*amp_factor)))) OK <- 1

                    ## If we have a really big amp_factor, it means we're not getting anywhere and
                    ## sigma2fshat is probably tending to zero.
                    if(amp_factor > 1e9) {
                        OK <- 1
                    }
                }

                if(amp_factor > 1e9) {
                    sigma2fs_new <- 0  # fix sigma2fshat to zero since we couldn't estimate it
                } else {
                    ## Otherwise find the root of the equation with the sought initial conditions
                    sigma2fs_new <- stats::uniroot(f = J,
                                                   interval = c(sigma2fs/amp_factor,
                                                                sigma2fs*amp_factor))$root
                }
            }
            
            D <- sigma2fs_new*Sm@Vfs + Sm@Ve  # total data variance-covariance
            if(isDiagonal(D)) {               # inverse of D (as above)
                D <- Diagonal(x=D@x)          # cast to Diagonal
                Dinv <- solve(D)
            } else {
                Dinv <- chol2inv(chol(D))
            }
            alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*%      # alpha GLS estimate (as above)
                t(Sm@X) %*% Dinv %*% (Sm@Z - Sm@S %*% mu_eta)
        } else {
            ## Here we are in the homoscedastic (diagonal) case and we can
            ## solve for alpha independently of sigma2fshat
            alpha <- solve(t(Sm@X) %*% Sm@X) %*%      # alpha GLS estimate
                t(Sm@X) %*% (Sm@Z - Sm@S %*% mu_eta)
            resid <- Sm@Z - Sm@X %*% alpha            # residual
            if(est_sigma2fs) {
                Omega_diag <- Omega_diag1 -               # just compute Omega once
                    2*diag2(Sm@S %*% mu_eta, t(resid)) +
                    diag2(resid,t(resid))
                Omega_diag <- Diagonal(x=Omega_diag)

                ## Closed-form solution for sigma2fs (see vignette)
                sigma2fs_new <- 1/b[1]*(sum(Omega_diag)/length(Sm@Z) - a[1])
                if(sigma2fs_new < 0) {  # If we get less than zero because of numeric instability
                    sigma2fs_new = 0    # just fix to zero
                }
            }            
        }
    }


    ## If we do NOT have any fine-scale variation (e.g., estimated to zero in previous iteration)
    if(all(diag(Sm@Vfs) == 0)) {
        alpha <- solve(t(Sm@X) %*% solve(Sm@Ve) %*% Sm@X) %*% t(Sm@X) %*%  # just find GLS
            solve(Sm@Ve) %*% (Sm@Z - Sm@S %*% mu_eta)
        if(est_sigma2fs) sigma2fs_new <- 0                                 # and keep sigma2fs at zero
    }

    ## Update SRE model with estimated quantities
    Sm@Khat <- K
    Sm@Khat_inv <- Khat_inv
    Sm@alphahat <- alpha
    Sm@sigma2fshat <- sigma2fs_new

    ## Return SRE model
    return(Sm)
}

## This routine updates the covariance matrix of the random effects
.update_K <- function(Sm,method="unstructured",
                      S_eta= NULL,mu_eta = NULL,
                      lambda = 0) {

    if (is.null(S_eta)) S_eta <- Sm@S_eta      # Conditional covariance matrix of random effects
    if (is.null(mu_eta)) mu_eta <- Sm@mu_eta   # Conditional mean of random effects

    if(method == "unstructured") {
        ## If K is unstructured, then the update is trivial, see vignette Section 2.2
        ## I allow for some regularisation through lambda should this be deemed required
        ## (This is useful for when we have lots of basis and few data points)
        K <- .regularise_K(Sm, lambda = lambda)
    } else if (method == "block-exponential") {
        ## If K is block exponential (blocked by resolution) then
        ## we need to find the (i) precision, (ii) spatial length scale, and
        ## (iii) temporal length scale by resolution
        all_res <- count_res(Sm)               # number of resolutions
        eta2 <- lapply(1:nrow(all_res),function(i) {
            ## find which indices correspond to these basis functions
            idx <- which(data.frame(Sm@basis)$res == i)
            S_eta[idx,idx] +
                tcrossprod(mu_eta[idx])
        })

        ## (i) Find the precision associated with each resolution
        omega <- lapply(1:nrow(all_res),       # for each resolution
                        function(i) {
                            ## number of basis functions in i-th resolution
                            ni <- all_res[i,]$n

                            ## find which indices correspond to these basis functions
                            idx <- which(data.frame(Sm@basis)$res == i)

                            ## # find the current CORRELATION matrix associated with this resolution
                            Ki <- Sm@Khat[idx,idx]/Sm@Khat[idx[1],idx[1]]

                            ## Compute INVERSE CORRELATION matrix associated with this resolution
                            Ki_inv <- chol2inv(chol(Ki))

                            ## The precision is given by n / tr(Kinv %*% (S_eta + mu.mu'))
                            ni / sum(diag2(Ki_inv,eta2[[i]]))
                        })

        ## (ii,iii) Likelihood function for spatial/temporal length scales
        f_tau <- function(tau_i,i) {   # tau_i are the scales, i is the resolution
            if(any(tau_i <= 1e-10)) {  # do not let any of the taus be too small
                Inf
            } else {

                ## Find which bases are at this resolution
                idx <- which(data.frame(Sm@basis)$res == i)

                ## Since we're block exponential, the correlation matrix is simply
                ## computed from the distances using the appropriate decay parameters
                if(is(Sm@basis,"TensorP_Basis")) {
                    ## If we have a tensor basis then construct Ki using the Kronecker product
                    Ki1 <- exp(-Sm@D_basis$Basis2[[1]]/tau_i[2])  # temporal part
                    Ki2 <- exp(-Sm@D_basis$Basis1[[i]]/tau_i[1])  # spatial part
                    ## time runs slowest (and only one time resolution),  space runs fastest
                    Ki <- kronecker(Ki1,Ki2)

                    ## Compute the inverse correlation matrix
                    Qi1 <- chol2inv(chol(Ki1))
                    Qi2 <- chol2inv(chol(Ki2))
                    Ki_inv <- kronecker(Qi1,Qi2)

                    ## Compute log determinant
                    R1 <- chol(Qi1)
                    R2 <- chol(Qi2)
                    det_part <- 0.5*(nrow(R2)*logdet(R1) + nrow(R1)*logdet(R2))

                    ## Compute the log=likelihood. There doesn't seem to be a way to
                    ## simplify this using the Kronecker product
                    -as.numeric(det_part - omega[[i]]/2*sum(diag2(Ki_inv,eta2[[i]],symm=TRUE)))

                } else {
                    ## Just spatial, from distances between centroid
                    Ki <- exp(-Sm@D_basis[[i]]/tau_i)

                    ## Compute the inverse correlation matrix
                    Ki_inv <- chol2inv(chol(Ki))

                    ## Compute the log-likelihood
                    -as.numeric(0.5*determinant(Ki_inv)$modulus -
                                    omega[[i]]/2*sum(diag2(Ki_inv,eta2[[i]],symm=TRUE)))
                }

            }
        }

        ## (ii,iii) GRADIENT of the likelihood function for spatial/temporal length scales
        gr_f_tau <- function(tau_i,i) {
            idx <- which(Sm@basis@df$res == i)  # tau_i are the scales, i is the resolution

            if(is(Sm@basis,"TensorP_Basis")) {
                Ki1 <- exp(-Sm@D_basis$Basis2[[1]]/tau_i[2])  # temporal part
                Ki2 <- exp(-Sm@D_basis$Basis1[[i]]/tau_i[1])  # spatial part
                Ki <- kronecker(Ki1,Ki2)                      # Kronecker of the two

                ## Compute the inverse correlation matrix
                Qi1 <- chol2inv(chol(Ki1))
                Qi2 <- chol2inv(chol(Ki2))
                Ki_inv <- kronecker(Qi1,Qi2)

                ## d(X kron Y) = dX kron Y + X cron dY. Compute these below
                dKi <- kronecker(Ki1,(Sm@D_basis$Basis1[[i]]/(tau_i[1]^2))*Ki2)
                dKit <- kronecker((Sm@D_basis$Basis2[[1]]/(tau_i[2]^2))*Ki1,Ki2)

            } else {
                ## If only spatial then just compute derivative of exponential
                Ki <- exp(-Sm@D_basis[[i]]/tau_i)
                dKi <- (Sm@D_basis[[i]]/(tau_i^2))*exp(-Sm@D_basis[[i]]/tau_i)
                Ki_inv <- chol2inv(chol(Ki))  # inverse
            }

            ## derivative of log-likelihood w.r.t tau_1 (spatial)
            tau_i1 <- -(-0.5*sum(diag2(dKi,Ki_inv)) +
                            0.5*omega[[i]]*sum(diag2(eta2[[i]]%*% Ki_inv,
                                                     dKi %*% Ki_inv)))
            tau_i1 <- as.numeric(tau_i1)

            if(length(tau_i) == 1) {  # Then we just have space
                return(tau_i1)
            } else {                  # We have time aswell

                ## derivative of log-likelihood w.r.t tau_2 (temporal)
                tau_i2 <-  -(-0.5*sum(diag2(dKit,Ki_inv)) +
                                 0.5*omega[[i]]*sum(diag2(eta2[[i]]%*% Ki_inv,
                                                          dKit %*% Ki_inv)))
                tau_i2 <- as.numeric(tau_i2)

                ## Return both derivatives
                return(c(tau_i1,tau_i2))
            }
        }

        ## Find the maximum spatial distance between centroids of all basis functions. This is used for initialisation
        max_l <- max(unlist(Sm@D_basis[[1]]))

        ## Below we actually estimate the parameters
        ## For each resolution
        tau <- lapply(1:nrow(all_res),
                      function(i) {
                          ## Find the basis functions for this resolution
                          idx <- which(Sm@basis@df$res == i)

                          ## Compute the correlation matrix
                          Ki <- Sm@Khat[idx,idx]/Sm@Khat[idx[1],idx[1]]

                          ## If we are in space-time
                          if(is(Sm@basis,"TensorP_Basis")) {
                              ## Extract previous estimate from current covariance matrix.
                              ## If zero (e.g., initial matrix is the identity), then pin to 1e-9
                              ## If only one basis function at this resolution then don't
                              ## attempt to estimate
                              par_init <- ifelse(all_res$n[i]>1,
                                                 max(-Sm@D_basis$Basis1[[i]][1,2]/log(Ki[1,2]),1e-9),
                                                 1e-9)


                              ## Same as above but for temporal (assume we have always more than one temporal basis function)
                              par_init[2] <- max(-Sm@D_basis$Basis2[[1]][1,2]/log(Ki[1,1+count_res(Sm@basis@Basis1)$n[i]]),1e-9) ## time

                              ## If we clamped the temporal length scale then set it initially to 1
                              if(par_init[2] == 1e-9) par_init[2] <- 1
                          } else {
                              ## As above but just for space
                              par_init <- ifelse(all_res$n[i]>1,
                                                 max(-Sm@D_basis[[i]][1,2]/log(Ki[1,2]),1e-9),
                                                 1e-9)
                          }

                          ## If we clamped the spatial length scale then set it initially to max(length) / 10
                          if(par_init[1] == 1e-9) par_init[1] <- max_l/10

                          ## Suppress warnings in case we hit max-iterations. If it hasn't converged we would be
                          ## in a GEM settings which is still OK
                          suppressWarnings(optim(par = par_init,
                                                 fn = f_tau,
                                                 gr = gr_f_tau,
                                                 i=i,control=list(maxit=100L,reltol=1e-4))$par)
                      })

        ## Reconstruct the K matrix based on above parameter estimates
        K <- lapply(1:nrow(all_res),
                    function(i) {
                        if(is(Sm@basis,"TensorP_Basis")) {
                            Ki <- kronecker(exp(-Sm@D_basis$Basis2[[1]]/tau[[i]][2]),
                                            exp(-Sm@D_basis$Basis1[[i]]/tau[[i]][1]))/omega[[i]]
                        } else {
                            Ki <- exp(-Sm@D_basis[[i]]/tau[[i]])/omega[[i]]
                        }

                    })

        ## Since we are in block diagonal mode we can just block-diagonalise across resolutions
        K <- do.call("bdiag",K)

        ## Now, if we have space AND time, block diagonalising by resolution is not correct
        ## as we have the following indices (res1t1....res1tN,res2t1,...,res2tN,...)
        ## This can be corrected by seeing how the indices were in the original data frame
        ## (which were correct by construction), and then permuting the K matrix using
        ## and internal function reverse_permute
        idx_all <- unlist(lapply(1:nrow(all_res),
                                 function(i) which(Sm@basis@df$res == i)))

        # reverse_permute rearranges the order of time/resolution when we have tensor products
        # When we don't have tensor product idx_all and 1:nrow(K) are the same so nothing changes
        K <- reverse_permute(K,idx_all)

        ## If user wants verbose output show estimates
        if( opts_FRK$get("verbose") > 0) {
            cat("  Estimates of omega: ",unlist(omega),"  ")
            cat("  Estimates of tau: ",unlist(tau),"  ")
        }

    }
    
    ## Return the estimated matrix
    return(K)
}

## The function below regularises the K matrix when the K_type is "unstructured"
.regularise_K <- function(Sm,S_eta= NULL,mu_eta = NULL, lambda = 0) {

    if (is.null(S_eta)) S_eta <- Sm@S_eta      # extract from SRE model if not supplied
    if (is.null(mu_eta)) mu_eta <- Sm@mu_eta   # extract from SRE model if not supplied

    if(any(lambda > 0)) {  # if at least one lambda > 0

        ## If we have just one regulatisation parameter for all resolutions
        if(length(lambda) == 1) {
            reg_matrix <- lambda*Diagonal(nrow(S_eta)) # reg. matrix = lambda*I
        } else {
            ## If we have one regularisation parameter per resolution then the reg. matrix
            ## is diagonal but not proportional to the identity matrix
            ## We use the data frame returned by count_res which has the resolution number
            ## in the first column and the number of basis in the second column
            reg_matrix <- Diagonal(x = do.call("c",
                                               apply(count_res(Sm),1,
                                                     function(x) rep(lambda[x[1]],x[2]))))
        }

        ## Update K but this time regularising
        Q <- chol2inv(chol(S_eta + tcrossprod(mu_eta))) + reg_matrix
        K <- chol2inv(chol(Q))
    } else {
        ## If there is no regularisation then use the following simple update (see vignette for details)
        K <- S_eta + tcrossprod(mu_eta)
    }

    ## Return K
    K
}


# ---- TMB fitting functions ----


## Fitting stage of non-Gaussian FRK (more generally, for method  = 'TMB').
.TMB_fit <- function(object, optimiser, known_sigma2fs, taper, ...) {
  
  ## Some matrices evaluated at observed BAUs only: 
  C_O <- .constructC_O(object) 
  X_O <- .constructX_O(object) 
  S_O <- .constructS_O(object) 
  
  info_fit <- list()      
  info_fit$time <- system.time({
    
  info_fit$method <- "TMB" 
  
  ## If we are using a precision matrix formulation, determine if we can use the
  ## neighbour formulation, or if we need to use the precision-block-exponential.
  ## This depends on whether the *spatial* basis is regular, or not.
  K_type <- object@K_type

  if (K_type == "precision") {
    sp_basis <- if (is(object@basis,"TensorP_Basis")) object@basis@Basis1 else object@basis
    if (!sp_basis@regular || is(manifold(sp_basis), "sphere")) {
      K_type <- "precision-block-exponential"
    } else {
      K_type <- "neighbour"
    }
  }
  
  if (is.null(taper) && (K_type %in% c("block-exponential", "precision-block-exponential"))) {
    cat("The argument taper was not specified. Since we are using method = 'TMB' 
         with either i) a covariance matrix (K_type = 'block-exponential') 
         or ii) irregular basis functions (object@basis@regular = 0) or iii) a 
         non-plane() manifold, we must use tapering for computational reasons. 
         Setting taper = 3.\n")
    taper <- 3
    info_fit$taper <- taper
  }
  
  ## Parameter and data preparation for TMB
  parameters <- .TMB_initialise(object, K_type = K_type, C_O = C_O, X_O = X_O, S_O = S_O)
  data <- .TMB_data_prep(object, sigma2fs_hat = exp(parameters$logsigma2fs), 
                         K_type = K_type, taper = taper, 
                         C_O = C_O, X_O = X_O, S_O = S_O)
  
  ## If we are estimating a unique fine-scale variance at each spatial BAU, 
  ## simply replicate sigma2fs ns times. 
  ns <- dim(object@BAUs)[1]
  if (object@fs_by_spatial_BAU) {
    data$sigma2fs_hat <- rep(data$sigma2fs_hat, ns)
    parameters$logsigma2fs <- rep(parameters$logsigma2fs, ns)
  }
  
  ## Fix sigma2fs to the known value provided by the user (if provided). 
  if (!is.null(known_sigma2fs)) {
    data$fix_sigma2fs <- 1
    data$sigma2fs_hat <- known_sigma2fs
    parameters$logsigma2fs <- log(known_sigma2fs) 
  }
  
  ## Don't want to pass in variance components that are "too small" or "too big". 
  ## This has caused TMB to explode and cause R to crash in the past, with no 
  ## explanation as to why the crash occured.
  parameters$logsigma2 <- pmin(pmax(parameters$logsigma2, -4), 8)
  parameters$logtau <- pmin(pmax(parameters$logtau, -4), 8)
  parameters$logdelta <- pmin(pmax(parameters$logdelta, -4), 8)
  parameters$logsigma2_t <- pmin(pmax(parameters$logsigma2_t, -4), 8)
  
  ## Checks to catch catastrophic errors. 
  ## If we allow nonsensical values into TMB, R may crash without providing
  ## an informative warning. These checks will hopefully ensure that will not happen. 
  if (any(sapply(data, function(x) any(length(x) == 0) || any(is.na(x)) || any(is.null(x)))))
    stop("Something has gone wrong in the data preparation for TMB: Some entries are numeric(0), NA, or NULL. Please contact the package maintainer.")
  if (any(sapply(parameters, function(x) any(length(x) == 0) || any(is.na(x)) || any(is.null(x)))))
    stop("Something has gone wrong in the parameter initialisation for TMB: Some entries are numeric(0), NA, or NULL. Please contact the package maintainer.")
  if (any(data$nnz < 0) || any(data$col_indices < 0) || any(data$row_indices < 0))
    stop("Something has gone wrong in construction of the precision matrix of the basis-function coefficients: We have negative row-indices, col-indices, or total non-zeros: Please contact the package maintainer. ")
  if (!all.equal(length(data$x), length(data$col_indices), length(data$row_indices), sum(data$nnz))) 
    stop("Something has gone wrong in construction of the precision matrix of the basis-function coefficients: The number of row-indices, col-indices, or non-zeros is inconsistent. Please contact the package maintainer. ")
  if(!all.equal(length(object@Z), length(data$Z), nrow(data$C_O), nrow(data$X_O) , nrow(data$S_O)))
    stop("Something has gone wrong in the data preparation for TMB: The dimensions of the C, X, or S matrix is inconsistent with the number of observations. Please contact the package maintainer.")
  if(!all.equal(nbasis(object@basis), max(data$row_indices + 1), max(data$col_indices + 1), sum(data$r_si)))
    stop("Something has gone wrong in the data preparation for TMB: The number of basis functions and the matrix indices are inconsistent. Please contact the package maintainer.")
  if (!(K_type %in% c("neighbour", "block-exponential", "precision-block-exponential")))
    stop("Internal error: K_type is not one of neighbour, block-exponential, or precision-block-exponential. Please contact the package maintainer.")


  ## TMB model compilation
  obj <- MakeADFun(data = data,
                   parameters = parameters,
                   random = c("random_effects"),
                   DLL = "FRK", 
                   silent = !opts_FRK$get("verbose")) # hide the gradient information during fitting
  
  ## The following means we want to print every parameter passed to obj$fn.
  obj$env$tracepar <- opts_FRK$get("verbose")

  
  cat("Optimising with TMB...\n")

  ## The optimiser should have arguments: start, objective, gradient. 
  ## The remaining arguments can be whatever.
  fit <- optimiser(obj$par, obj$fn, obj$gr, ...)
  
  info_fit$iterations  <- fit$iterations
  info_fit$convergence <- fit$convergence
  info_fit$message     <- fit$message
  
  cat("Optimisation completed.\n")

  
  # ---- Joint precision/covariance matrix of random effects ----
  
  cat("Extracting estimates of the parameters and the joint precision matrix of the random effects from TMB...\n")
  
  ## Extract parameter and random effect estimates
  par <- obj$env$last.par.best
  estimates <- split(par, names(par)) # convert to named list object
  
  ## TMB treats all parameters (fixed effects, variance components, 
  ## and random effects) as random quantities, and so the joint precision 
  ## matrix obtained using sdreport(obj, getJointPrecision = TRUE) contains 
  ## the precision matrix for fixed and random effects. 
  ## However, we assume the regression parameters and variance components are 
  ## fixed effects, NOT random quantities, and so they should not have a 
  ## randomness associated to them. We overcome this by considering the fixed 
  ## effects and parameters as random during the fitting process, and then 
  ## post-fitting we condition on theta = \hat{theta}, the ML estimate.
  ## By conditioning on \hat{theta}, we can consider the precision matrix of  
  ## the random effects in isolation.
  s <- length(estimates$random_effects) # Number of random effects
  
  if (object@simple_kriging_fixed) {
    ## Extract the precision matrix of the random effects only
    Q_posterior <- obj$env$spHess(par = obj$env$last.par.best, random = TRUE)
  } else {
    ## if we wish to extract the uncertainty of the parameters and fixed effects
    ## in order to do universal kriging (or, at least give us the OPTION to do 
    ## it during the prediction stage), We need to use sdreport().
    ## This can be computationally demanding, which is why we provide the
    ## argument simple_kriging_fixed.
    Q_posterior <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
    
    ## We will only retain the uncertainty in the fixed effects
    ## (i.e., in alpha), and not the parameters.
    retain_idx  <- rownames(Q_posterior) %in% c("alpha", "random_effects") 
    Q_posterior <- Q_posterior[retain_idx, retain_idx]
  }
  
  cat("Extraction completed.\n")
  
  ## Update the slots of object
  ## Convert to Matrix as these SRE slots require class "Matrix"
  r  <- nbasis(object)
  mstar <- ncol(C_O)
  object@alphahat <- as(estimates$alpha, "Matrix")
  object@mu_eta   <- as(estimates$random_effects[1:r], "Matrix")
  if (object@include_fs) {
    object@mu_xi  <- as(estimates$random_effects[(r+1):(r + mstar)], "Matrix")
  } else {
    object@mu_xi  <- as(rep(0, mstar), "Matrix")
  }
  
  object@sigma2fshat <- unname(exp(estimates$logsigma2fs))
  object@Q_posterior <- Q_posterior
  object@phi <- unname(exp(estimates$logphi))
  
  ## Log-likeihood (negative of the negative-log-likelihood)
  object@log_likelihood <- -obj$fn() # could also use -fit$objective
  
  ## If zero fine-scale variation detected just make sure user knows.
  ## This can be symptomatic of poor fitting
  if(object@sigma2fshat == 0) {
    info_fit$sigma2fshat_equal_0 <- 1
    if(opts_FRK$get("verbose") > 0)
      message("sigma2fs is being estimated to zero.
      This might because of an incorrect binning procedure or because too much 
      measurement error is being assumed (or because the latent field is indeed 
              that smooth, but unlikely).")
  } else {
    info_fit$sigma2fshat_equal_0 <- 0
  }
  
  })
  
  object@info_fit <- info_fit

  return(object)
}



## Initialise the fixed effects, random effects, and parameters for method = 'TMB'
.TMB_initialise <- function(object, K_type, C_O, X_O, S_O) {   
  
  cat("Initialising parameters and random effects...\n")
  
  nres    <- max(object@basis@df$res) 
  r       <- object@basis@n  
  
  # ---- Estimate Y over the observed BAUs ----
  
  Y_O <- .compute_Y_O(object, C_O = C_O)
  
  # ---- Parameter and random effect initialisations ----
  
  ## Now that we have a crude estimate of Y_O, we may initialise the variance
  ## components and random effect
  
  l <- list() # list of initial values
  
  ## i. Fixed effects alpha (OLS solution)
  l$alpha <- solve(t(X_O) %*% X_O) %*% t(X_O) %*% Y_O # OLS solution
  ## ii. Variance components
  ## Dispersion parameter depends on response; some require it to be 1. 
  if (object@response %in% c("poisson", "binomial", "negative-binomial")) {
    l$phi <- 1
  } else if (object@response == "gaussian") {
    l$phi <- mean(diag(object@Ve))
  } else {
    ## Use the variance of the data as our estimate of the dispersion parameter.
    ## This will almost certainly be an overestimate, as the mean-variance 
    ## relationship is not considered.
    l$phi <- var(as.numeric(object@Z))
  }

  
  
  ## Variance components of the spatial basis-function coefficients
  l$sigma2  <- var(as.vector(Y_O)) * (0.1)^(0:(nres - 1))
  l$tau     <- (1 / 3)^(1:nres)
  if (K_type != "block-exponential") {
    # l$sigma2   <- 1 / exp(l$sigma2) 
    # l$tau      <- 1 / exp(l$tau)   
    l$sigma2   <- 1 / l$sigma2
    l$tau      <- 1 / l$tau
  } 
  
  
  if (K_type == "separable") {
    ## Separability means we have twice as many spatial basis function variance
    ## components. So, just replicate the already defined parameters.
    l$sigma2 <- rep(l$sigma2, 2)
    l$tau    <- rep(l$tau, 2)
    l$logdelta <- 1 # Dummy value
  } else if (K_type == "precision-block-exponential") {
    ## Precision exp requires one extra parameter per resolution
    l$logdelta <- rnorm(nres)
  } else {
    l$logdelta <- 1 # Dummy value
  }
  
  ## Variance components of temporal basis-function coefficients
  l$sigma2_t    <- 1
  l$rho_t       <- 0.1
  
  for (iteration_dummy in 1:5) {
    
    ## iii. Basis function random-effects 
    regularising_weight <- if (!is.null(l$sigma2fs)) l$sigma2fs else l$sigma2[1] 
    
    QInit <- .sparse_Q_block_diag(object@basis@df, 
                                  # kappa = exp(l$sigma2), 
                                  # rho = exp(l$tau))$Q      
                                  kappa = l$sigma2, 
                                  rho = l$tau)$Q
    
    ## Matrix we need to invert
    mat <- Matrix::t(S_O) %*% S_O / regularising_weight + QInit 
    
    ## Invert the matrix
    mat_inv <- tryCatch(expr = {
      if (r > 4000) { # Avoid full inverse if we have too many basis functions
      mat_L <- sparseinv::cholPermute(Q = mat) # Permuted Cholesky factor
      ## Sparse-inverse-subset (a proxy for the inverse)
      sparseinv::Takahashi_Davis(Q = mat, cholQp = mat_L$Qpermchol, P = mat_L$P)
    } else {
      solve(mat) 
    }}, 
    error = function(w){
      cat("Initialisation phase: could not invert mat, just using diag(r)\n")
      diag(r)
      }
    )

    if (regularising_weight == 0) { # avoid division by zero
      warning("In initialisation stage, the regularising_weight is 0; setting it to 1. This is probably not an issue, but feel free to contact the package maintainer.")
      regularising_weight <- 1
    }
    
    ## MAP estimate of eta
    l$eta  <- (1 / regularising_weight) * mat_inv %*% Matrix::t(S_O)  %*% (Y_O - X_O %*% l$alpha)
    
    ## iv. Observed fine-scale random effects xi_O
    l$xi_O <- Y_O - X_O %*% l$alpha - S_O %*% l$eta
    
    ## v. Fine-scale variance
    l$sigma2fs <- var(as.vector(l$xi_O)) 
  }
  
  ## Return list of parameter initialisations for TMB
  transform_minus_one_to_one_inverse <- function(x) -0.5 * log(2 / (x + 1) - 1)
  return(list(
    alpha = as.vector(l$alpha),
    logphi = log(l$phi),
    logsigma2 = log(l$sigma2),
    logtau = log(l$tau),
    logdelta = l$logdelta,
    logsigma2_t = log(l$sigma2_t),  
    frho_t = transform_minus_one_to_one_inverse(l$rho_t),
    logsigma2fs = log(l$sigma2fs),
    random_effects = c(as.vector(l$eta), if(object@include_fs) as.vector(l$xi_O))
  ))
}


.compute_Y_O <- function(object, C_O) {
  
  Z       <- as.vector(object@Z)         
  k_Z     <- object@k_Z       
  k_BAU_O <- object@k_BAU_O 
  
  ## Create the relevant link functions. When a probability parameter is 
  ## present in a model and the link-function is appropriate for modelling 
  ## probabilities (i.e., the link function maps to [0, 1]), we may use 
  ## hierarchical linking to first link the probability parameter to the 
  ## Gaussian Y-scale, and then the probability parameter to the conditional 
  ## mean at the data scale. In other situations, we simply map from Y to the mean.
  if (object@response %in% c("binomial", "negative-binomial") & object@link %in% c("logit", "probit", "cloglog")) {
    f     <- .link_fn(kind = "prob_to_Y", link = object@link)
    h     <- .link_fn(kind = "mu_to_prob", response = object@response)
  } else {
    g     <- .link_fn(kind = "mu_to_Y", link = object@link) 
  }
  
  ## Create altered data to avoid the problems of applying g() to Z.
  ## This altered data is used only during the initialisation stage.
  Z0 <- Z
  if (object@link %in% c("log", "sqrt")) {
    Z0[Z <= 0] <- 0.1      
  } else if (object@response == "negative-binomial" & object@link %in% c("logit", "probit", "cloglog")) {
    Z0[Z == 0]   <- 0.1
  } else if (object@response == "binomial" & object@link %in% c("logit", "probit", "cloglog")) {
    Z0 <- Z + 0.1 * (Z == 0) - 0.1 * (Z == k_Z)
  } else if (object@link %in% c("inverse-squared", "inverse")) {
    Z0[Z == 0] <- 0.05
  } 
  
  # ---- Estimate mu_Z, mu_O, and Y_O ----
  
  ## First, we estime mu_Z with the (adjusted) data
  mu_Z <- Z0
  

  ## Use mu_Z to construct mu_O
  mu_O <- .compute_mu_O(object, C_O, mu_Z)
  
  ## For some link functions, mu_0 = 0 causes NaNs; set these to a small positive value.
  ## The size parameter being 0 also causes NaNs.
  k_BAU_O[k_BAU_O == 0] <- 1
  mu_O <- mu_O + 0.05 * (mu_O == 0) - 0.05 * (mu_O == k_BAU_O)
  
  
  ## Transformed data: convert from data scale to Gaussian Y-scale.
  if (object@response %in% c("binomial", "negative-binomial") & object@link %in% c("logit", "probit", "cloglog")) {
    Y_O <- f(h(mu_O, k_BAU_O)) 
  } else if (object@response == "negative-binomial" & object@link %in% c("log", "sqrt")) {
    Y_O <- g(mu_O / k_BAU_O) 
  } else {
    Y_O <- g(mu_O)
  } 
  
  return(Y_O)
}


.compute_mu_O <- function(object, C_O, mu_Z) {

  m       <- nrow(C_O)
  mstar   <- ncol(C_O)
  obs_BAUs_df <- object@BAUs@data[object@obsidx, ]

  mu_O <- vector(mode = "list", length = mstar)
  for (Bj in 1:m) {        # for each observation (obs.) j

    w_j <- C_O[Bj, ]              # extract the incidence matrix weights for obs. j
    idx <- which(w_j > 0)         # find the BAU indices associated with obs. j
    
    ## If the data is bin. or neg.-bin., use the BAU level size parameters to 
    ## interpolate mu_O from mu_Z. (Note that this step is the reason why we enforce
    ## the weights of the incidence matrix to be 1 if the data is bin. or neg.-bin.)
    ## Otherwise, just use the elements (weights) of C_Z to interpolate mu_O from mu_Z.
    if (object@response %in% c("binomial", "negative-binomial")) {
      w_j <- obs_BAUs_df$k_BAU[idx]
    } else {
      w_j <- w_j[idx]
    }

    ## Interpolate mu_Z[j] to each BAU associated with obs. j
    mu_O_j <- w_j / sum(w_j) * mu_Z[Bj] 
    for (i in 1:length(idx)) {  # for each BAU associated with obs. j,
      mu_O[[idx[i]]] <- c(mu_O[[idx[i]]], mu_O_j[i]) # assign mu_O_j[i] to BAU i
    }
  }
  
  ## If we have overlapping data supports, some entries in the list mu_O will 
  ## have a length greater than 1. If the data is bin. or neg.-bin., take the 
  ## minimum for each entry, because we must have mu_O[i] <= k_BAU_O[i] for all 
  ## i = 1...N. Otherwise, just take the average. 
  if (object@response %in% c("binomial", "negative-binomial")) {
    mu_O <- sapply(mu_O, min)
  } else {
    mu_O <- sapply(mu_O, mean)
  }
  
  return(mu_O)
}

.TMB_data_prep <- function (object, sigma2fs_hat, K_type, taper, C_O, X_O, S_O) {
  
  cat("Preparing data for TMB...\n")
  
  obsidx <- observed_BAUs(object)       # Indices of observed BAUs
  
  ## measurement error variance (NB: this is a vector)
  sigma2e <- if (object@response == "gaussian") diag(object@Ve) else -1
  
  ## Data passed to TMB which is common to all
  data <- list(Z = as.vector(object@Z),  # Binned data
               X_O = X_O, S_O = S_O, C_O = C_O,
               K_type = K_type, response = object@response, link = object@link,
               k_BAU_O = object@k_BAU_O, k_Z = object@k_Z,         
               temporal = as.integer(is(object@basis,"TensorP_Basis")), 
               fs_by_spatial_BAU = object@fs_by_spatial_BAU, sigma2e = sigma2e, 
               BAUs_fs = object@BAUs$fs[obsidx])
  
  ## Define the number of temporal basis function (r_t), and number of spatial BAUs (ns).
  ns <- dim(object@BAUs)[1]
  if (data$temporal) {
    spatial_dist_matrix <- object@D_basis[[1]]
    spatial_basis <- object@basis@Basis1  
    data$r_t <- object@basis@Basis2@n
  } else {
    spatial_dist_matrix <- object@D_basis
    spatial_basis <- object@basis 
    data$r_t <- 1
  }
  
  data$spatial_BAU_id <-  (obsidx - 1) %% ns
  data$r_si <- as.vector(table(spatial_basis@df$res))
  
  ## Data which depend on K_type: provide dummy data (can't provide nothing)
  ## and only change the ones we actually need.
  data$beta <- data$nnz <- data$row_indices <- data$col_indices <- 
    data$x <- data$n_r <- data$n_c  <- -1 
  
  if (K_type %in% c("block-exponential", "precision-block-exponential")) {
    tmp         <- .cov_tap(spatial_dist_matrix, taper = taper)
    data$beta   <- tmp$beta 
    R            <- as(tmp$D_tap, "dgTMatrix")
    data$nnz         <- tmp$nnz 
    data$row_indices <- R@i
    data$col_indices <- R@j
    data$x           <- R@x
    
  } else if (K_type == "neighbour") {
    tmp <- .sparse_Q_block_diag(spatial_basis@df, kappa = 0, rho = 1)
    R <- as(tmp$Q, "dgTMatrix")
    data$nnz         <- tmp$nnz
    data$row_indices <- R@i
    data$col_indices <- R@j
    data$x           <- R@x
    
  } else if (K_type == "separable") {
    ## Compute number of basis functions in each row (n_r) and each column (n_c)
    for (i in unique(spatial_basis@df$res)) {
      tmp <- spatial_basis@df[spatial_basis@df$res == i, ]
      data$n_r[i] <- length(unique(tmp$loc1))
      data$n_c[i] <- length(unique(tmp$loc2))
    }
  } 

  ## Create a data entry of sigma2fs_hat (one that will stay constant if we are 
  ## not estimating sigma2fs within TMB)
  data$sigma2fs_hat <- sigma2fs_hat
  ## fix sigma2fs during model fitting if all observations are associated with multiple
  ## BAUs, as TMB will explode if we don't. Otherwise, just estimate it with TMB.
  if (!any(tabulate(object@Cmat@i + 1) == 1)) {
    cat("There no observations are associated with a single BAU (i.e., all observations are associated with multiple BAUs). This makes the fine-scale variance parameter very difficult to estimate, so we will estimate it offline and fix for the remainder of model fitting; this estimate may be inaccurate.\n")
    data$fix_sigma2fs <- as.integer(1)
    if (object@fs_by_spatial_BAU) 
      stop("We do not allow each spatial BAU to have its own fine-scale variance parameter when there no observations associated with a single BAU (i.e., all observations are associated with multiple BAUs).")
  } else {
    data$fix_sigma2fs <- as.integer(0)
    if (!all(tabulate(object@Cmat@i + 1) == 1)) 
      cat("Some (but not all) observations are associated with multiple BAUs. Estimation of the fine-scale variance parameter will be done using TMB, but note that there should be a reasonable number of fine-unit observations so that TMB can get a handle of the fine-scale variance parameter.\n")
  } 
  
  data$include_fs   <- as.integer(object@include_fs)
  return(data)
}
