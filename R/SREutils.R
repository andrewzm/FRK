
SRE <- function(f,data,basis) {
    if(!is(f,"formula")) stop("f needs to be a formula.")
    if(!is(data,"Spatial")) stop("data needs to be of class Spatial (package sp).")
    if(!is(basis,"Basis")) stop("basis needs to be of class Basis (package FRK)")
    if(!("std" %in% names(data))) stop("data needs to contain a field 'std' denoting the observation error")
    if(!("fs" %in% names(data))) {
        warning("data should contain a field 'fs' containing a basis function for fine-scale variation. Setting basis function equal to one everywhere.")
        data$fs <- 1
    }
    if(!(all(data$fs > 0))) stop("fine-scale variation basis function needs to be nonnegative everywhere")

    L <- gstat:::gstat.formula(f,data=data)

    X <- as(L$X,"Matrix")
    Z <- Matrix(L$y)
    Ve <- Diagonal(x=data$std^2)
    Vfs <- Diagonal(x=data$fs)

    new("SRE",
        data=data,
        basis=basis,
        f = f,
        S = eval_basis(basis, s = coordinates(data)),
        Ve = Ve,
        Vfs = Vfs,
        X = X,
        Z = Z,
        mu_eta = Matrix(0,nbasis(basis),1),
        S_eta = Diagonal(x = rep(1,nbasis(basis))),
        alphahat = solve(t(X) %*% X) %*% t(X) %*% Z,
        Khat = Diagonal(n=nbasis(basis),x = var(Z[,1])),
        sigma2fshat = mean(diag(Ve)) / mean(diag(Vfs)))
}

SRE.fit <- function(SRE_model,n_EM = 100L, tol = 1e-5, method="EM",print_lik=F) {
    n <- nbasis(SRE_model)
    X <- SRE_model@X
    lk <- rep(0,n_EM)

    pb <- txtProgressBar(min = 0, max = n_EM, style = 3)
    for(i in 1:n_EM) {
        lk[i] <- loglik(SRE_model)
        SRE_model <- SRE.Estep(SRE_model)
        SRE_model <- SRE.Mstep(SRE_model)
        setTxtProgressBar(pb, i)
        if(i>2) if(lk[i] - lk[i-1] < tol) break
    }
    close(pb)

    if(i == n_EM) print("Maximum EM iterations reached")
    if(print_lik) plot(lk[1:i][-1])
    SRE_model
}

setMethod("SRE.predict",signature(Sm="SRE", pred_locs="Spatial",depname="character"),
          function(Sm,pred_locs,depname) {

              if(!("fs" %in% names(pred_locs))) {
                  warning("data should contain a field 'fs' containing a basis function for fine-scale variation. Setting basis function equal to one everywhere.")
                  pred_locs$fs <- 1
              }
              if(!(all(pred_locs$fs > 0))) stop("fine-scale variation basis function needs to be nonnegative everywhere")
              pred_locs[[depname]] <- median(Sm@data[[depname]])
              L <- gstat:::gstat.formula(Sm@f,data=pred_locs)
              X = as(L$X,"Matrix")
              S0 <- eval_basis(Sm@basis,pred_locs)
              pred_locs[[depname]] <- NULL

              alpha <- Sm@alphahat
              K <- Sm@Khat
              sigma2fs <- Sm@sigma2fshat
              D <- sigma2fs*Sm@Vfs + Sm@Ve
              Dinv <- solve(D)



              sig2_Vfs_pred <- Diagonal(x=sigma2fs*pred_locs$fs)
              if(is(pred_locs,"SpatialPointsDataFrame")) {

                    i_idx <- 1:nrow(Sm@Z)
                    j_idx <- which(pred_locs@data$id %in% Sm@data$id )
                    C <- sparseMatrix(i=i_idx,j=j_idx,x=1,dims=c(length(Sm@Z),nrow(pred_locs)))

                    LAMBDA <- bdiag(Sm@Khat,sig2_Vfs_pred)
                    LAMBDAinv <- chol2inv(chol(LAMBDA))
                    PI <- cBind(S0, .symDiagonal(n=nrow(pred_locs)))
                    Qx <- t(PI) %*% t(C) %*% solve(Sm@Ve) %*% C %*% PI + LAMBDAinv
                    temp <- cholPermute(Qx)
                    ybar <- t(PI) %*%t(C) %*% solve(Sm@Ve) %*% (Sm@Z - C %*% X %*% alpha)
                    x_mean <- cholsolve(Qx,ybar,perm=T,cholQp = temp$Qpermchol, P = temp$P)
                    Partial_Cov <- Takahashi_Davis(Qx,cholQp = temp$Qpermchol,P = temp$P)
                    x_margvar <- diag(Partial_Cov)

                    pred_locs[["mu"]] <- as.numeric(X %*% alpha + PI %*% x_mean)
                    pred_locs[["var"]] <- as.numeric(rowSums((PI %*% Partial_Cov)*PI))

#                     S_eta <- chol2inv(chol(crossprod(sqrt(Dinv) %*% Sm@S) + solve(K)))
#                     mu_eta <- S_eta %*% t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha)
#                     pred_locs[["mu"]] <- as.numeric(X %*% alpha + S0 %*% mu_eta)
#                     pred_locs[["var"]] <- as.numeric(rowSums((S0 %*% S_eta)*S0) +
#                                                          diag(sig2_Vfs_pred))
                } else if (is(pred_locs,"SpatialPolygonsDataFrame")) {
                    i_idx <- 1:nrow(Sm@Z)
                    j_idx <- over(Sm@data,as(pred_locs,"SpatialPolygons"))
                    C <- sparseMatrix(i=i_idx,j=j_idx,x=1,dims=c(length(Sm@Z),nrow(pred_locs)))

                    LAMBDA <- bdiag(Sm@Khat,sig2_Vfs_pred)
                    LAMBDAinv <- chol2inv(chol(LAMBDA))
                    PI <- cBind(S0, .symDiagonal(n=nrow(pred_locs)))
                    Qx <- t(PI) %*% t(C) %*% solve(Sm@Ve) %*% C %*% PI + LAMBDAinv
                    temp <- cholPermute(Qx)
                    ybar <- t(PI) %*%t(C) %*% solve(Sm@Ve) %*% (Sm@Z - C %*% X %*% alpha)
                    x_mean <- cholsolve(Qx,ybar,perm=T,cholQp = temp$Qpermchol, P = temp$P)
                    Partial_Cov <- Takahashi_Davis(Qx,cholQp = temp$Qpermchol,P = temp$P)
                    x_margvar <- diag(Partial_Cov)

                    pred_locs[["mu"]] <- as.numeric(x_mean[-(1:Sm@basis@n)])
                    pred_locs[["var"]] <- as.numeric(x_margvar[-(1:Sm@basis@n)])
                }
              pred_locs

          })




SRE.Estep <- function(Sm) {

    alpha <- Sm@alphahat
    K <- Sm@Khat
    sigma2fs <- Sm@sigma2fshat

    D <- sigma2fs*Sm@Vfs + Sm@Ve
    Dinv <- solve(D)

    S_eta <- chol2inv(chol(crossprod(sqrt(Dinv) %*% Sm@S) + solve(K)))
    mu_eta <- S_eta %*% t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha)

    Sm@mu_eta <- mu_eta
    Sm@S_eta <- S_eta

    Sm
}

SRE.Mstep <- function(Sm) {

    mu_eta <- Sm@mu_eta
    S_eta <- Sm@S_eta
    alpha_init <- Sm@alphahat
    sigma2fs_init <- Sm@sigma2fshat

    K <- S_eta + tcrossprod(mu_eta)
    alpha <- alpha_init
    sigma2fs <- sigma2fs_init
    converged <- FALSE



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
        Omega_diag <- diag2(Sm@S %*% S_eta,t(Sm@S)) +
            diag2(Sm@S %*% mu_eta %*% t(mu_eta), t(Sm@S)) -
            2*diag2(Sm@S %*% mu_eta, t(resid)) +
            diag2(resid,t(resid))
        Omega_diag <- Diagonal(x=Omega_diag)

        sigma2fs_new <- uniroot(f = J,interval = c(sigma2fs/100,sigma2fs*100))$root
        D <- sigma2fs_new*Sm@Vfs + Sm@Ve
        Dinv <- solve(D)
        alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*% t(Sm@X) %*% Dinv %*% (Sm@Z - Sm@S %*% mu_eta)

        if(max(sigma2fs_new / sigma2fs, sigma2fs / sigma2fs_new) < 1.001) converged <- TRUE
        sigma2fs <- sigma2fs_new
    }

    Sm@Khat <- K
    Sm@alphahat <- alpha
    Sm@sigma2fshat <- sigma2fs

    Sm
}



loglik <- function(Sm) {

    D <- Sm@sigma2fshat*Sm@Vfs + Sm@Ve
    Dinv <- solve(D)
    r <- Sm@Z - Sm@X %*% Sm@alphahat  - Sm@S %*% Sm@mu_eta

    as.numeric(-0.5*determinant(D)$modulus -
                   0.5 * t(r) %*% Dinv %*% r)
}

tr <- function(X) {
    sum(diag(X))
}

diag2 <- function(X,Y) {
    rowSums(X * t(Y))
}

