#' @title Get scale from practical range
#' 
#' @description Takes a practical range \code{l} and the smoothness of the Matern kernel \eqn{\nu} and returns the scale parameter \eqn{\kappa}.
#' @param l the practical range.
#' @param nu the smoothness parameter
#' @return numeric
#' @keywords Matern, practical range
#' @export
#' @seealso The inverse \code{\link{l_from_kappa}}
#' @examples
#' l = 100
#' nu = 1
#' kappa <- kappa_from_l(l,nu)
kappa_from_l <- function(l,nu) { return(sqrt(8*nu)/l)}

#' @title Get practical range from length scale
#' 
#' @description Takes a scale parameter \eqn{\kappa} and the smoothness of the Matern kernel \eqn{\nu} and returns the practical range \code{l}.
#' @param kappa the scale parameter
#' @param nu the smoothness parameter
#' @return numeric
#' @keywords Matern, practical range
#' @export
#' @seealso The inverse \code{\link{kappa_from_l}}
#' @examples
#' kappa = 0.02
#' nu = 1
#' l <- l_from_kappa(kappa,nu)
l_from_kappa <- function(kappa,nu) { return(sqrt(8*nu)/kappa)}


#' @title Precision matrix from SPDE
#' 
#' @description This function returns the precision matrix for a GMRF over a triangulation such that it is an approximation 
#' to a Gaussian field with a Matern kernel with pre-specified parameters.
#' @param M the mass matrix of the finite element basis
#' @param K the stiffness matrix of the finite element basis
#' @param nu the smoothness parameter of the Matern kernel
#' @param desired_prec a constant or a vector of length \code{n} where each element contains the 
#' precision (1 over the marginal standard deviation squared) at each vertex.
#' @param l a constant or a vector of length \code{n} where each element contains the 
#' practical range (distance at which correlation is 0.1) at each vertex. This should be spatially smooth. Informally, 
#' this should vary slower than the field itself, otherwise a good local approximation to the Matern representation is not guaranteed.
#' @return Object of class "dgCMatrix"
#' @keywords GMRF
#' @export
#' @examples
#' require(Matrix)
#' data(surf_fe)
#' mu <- matrix(0,nrow(surf_fe$p),1) 
#' Q <- Prec_from_SPDE_wrapper(surf_fe$M,
#'                             surf_fe$K,
#'                             nu = 2,
#'                             desired_prec = 1,
#'                             l = 100)
#' my_GMRF <- GMRF(mu=mu, Q=Q,name="my_first_GMRF")
Prec_from_SPDE_wrapper <- function(M,K,nu,desired_prec,l) {
  kappa_l = sqrt(8*nu)/l
  marg_prec <- marg_prec_from_kappa(kappa_l,nu)
  tau <- sqrt(desired_prec/marg_prec)
  Q <- Prec_from_SPDE(M,K,tau=tau,kappa=kappa_l,alpha=nu+1)
  return(Q)
}

#' @title The Matern function
#' 
#' @param r a vector or matrix of distances.
#' @param nu the smoothness parameter.
#' @param var the marginal variance
#' @param kappa the scale parameters
#' @return a vector or matrix after passing \code{r} through the Matern function.
#' @keywords Matern function, covariance function.
#' @export
#' @examples
#' x1 <- runif(10)
#' Sigma <- Matern(as.matrix(dist(x1)),nu=3/2,var=1,kappa=100)
Matern <- function(r=0:100,nu=3/2,var=1,kappa=0.1) {
  K <- var/((2^(nu-1))*gamma(nu))*(kappa*abs(r))^nu*besselK(kappa*abs(r),nu=nu)
  if (class(K) == "matrix")
    if(nrow(K) == ncol(K)) {
    diag(K) = var
  }
  return(K)
}

#' @title Return neighb list from precision structure
#' @description Takes a precision matrix and returns a list of vectors, where the \code{i}th vector contains \code{i} and the indices of the neighbours of the \code{i}th variable.
#' @param Q a sparse matrix
#' @return a list of vectors as described in the description
#' @export
#' @examples
#' G <- GMRF_RW()
#' nlist <- neighb_from_prec(getPrecision(G))
neighb_from_prec <- function(Q) {
  n <- dim(Q)[1]
  apply(matrix(1:n),1,function(x) { return(which(!(Q[x,]==0)))})
}


#' @title Return adjacency matrix from list of neighbours
#' @author Andrew Zammit Mangion
#' @description Creates a sparse adjacency matrix from list of vertex neighbours.
#' @param nei_list a list of vectors, where each vector contains the indices of the adjacent neighbours.
#' @return a sparse adjacency matrix (with zero on the diagonal)
#' @export
#' @examples
#' nei_list <- list(c(2,3),1,1)
#' adj.matrix <- adj_matrix_from_nei_list(nei_list)
adj_matrix_from_nei_list <- function(nei_list) {
  
  adj.matrix <- sparseMatrix(i= unlist(sapply(1:length(nei_list),function(i) {rep(i,length(nei_list[[i]]))})),
                             j = as.vector(unlist(nei_list)),
                             x=1)
  return(adj.matrix)
}


setMethod("sample_GMRF",signature="VAR_Gauss",function(G,L=NULL,reps=1,P=NULL) {
  G <- as(G,"GMRF")
  G@n <- nrow(G@Q)
  return(sample_GMRF(G,L,reps,P))
})

setMethod("sample_GMRF",signature="GMRF",function(G,L=NULL,reps=1,P=NULL) {
  n = G@n
  z <- matrix(rnorm(n*reps),n,reps)
  mu <- c(getMean(G))
  if(is.null(L)) {
    L <-t(chol(G@Q))
    P <-NULL
  }
  
  # Algorithm 2.4, Rue and Held
  if (is.null(P)) {
    v <- solve(t(L),z)
    x <- mu + v
  } else {
    v <-  P %*% solve(t(L),z)
    x <- mu + v
  }
  ifelse(reps==1, return(as.vector(x)), return(x))
})


### Loopy version
# setMethod("sample_GMRF",signature="GMRF",function(G,L=NULL,reps=1,P=NULL) {
#   n = G@n
#   x <- matrix(0,n,reps)
#   z <- matrix(rnorm(n*reps),n,reps)
#   mu <- getMean(G)
#   if(is.null(L)) {
#     L <-t(chol(G@Q))
#     P <-NULL
#   }
#   
#   # Algorithm 2.4, Rue and Held
#   if (is.null(P)) {
#     for (i in (1:reps)) {
#       v <- solve(t(L),z[,i])
#       x[,i] <- matrix(mu + v)
#     }
#   } else {
#     for (i in (1:reps)) {
#       v <- P %*% solve(t(L), z[,i])
#       x[,i] <- matrix(mu + v)
#     }
#   }
#   ifelse(reps==1, return(as.vector(x)), return(x))
# })

Prec_from_lattice <- function(Grid,ds) {
  n <- nrow(Grid)
  #### set up distance and neighbourhood (W, based on sharing a common border) matrices
  distance <-array(0, c(n,n))
  W <- as.matrix(dist(Grid))==ds
  n.neighbours <- as.numeric(apply(W, 1, sum))
  Q <- diag(n.neighbours) - W
  return(as(Q,"dgCMatrix"))
}


## Create a precision matrix from neighbourhood list
Prec_from_neighb <- function(neighb,intrinsic=1,precinc=1){
  
  num_v <- length(neighb)
  num_neighb <- lapply(neighb,length)
  
  if (intrinsic == 1)  {
    i_list <-  vector("list",num_v)    # row number index
    for (k in 1:num_v) {
      i_list[[k]] <- rep(k,num_neighb[[k]])
    }
    i <- unlist(i_list)
    j <- unlist(neighb)             # column number index
    
    z <- rep(-1, length(j))  # for intrinsic = 1 all non-zero off-diagonal elements are -1
    
    # Now add the diagonal elements (= number of neighbours
    i <- c(i,1:num_v)
    j <- c(j,1:num_v)
    zdiag <- unlist(num_neighb)
    #zdiag[zdiag == 0] = 1   # Make lonely elements conditionally independent with finite variance
    z <- c(z,zdiag)
  }
  
  if (intrinsic == 2)  {
    
    # Start off with the diagonal elements
    i1 <- 1:num_v
    j1 <- 1:num_v
    z1 <- rep(0,num_v)
    for (k in 1:num_v) {
      z1[k] <- num_neighb[[k]]^2 + num_neighb[[k]]
    }
    
    # Now add elements corresponding to the neighbours. Initialise, prob max 10 neighbours per node
    count <- 1
    i2 <- rep(0,num_v*10)
    j2 <- rep(0,num_v*10)
    z2 <- rep(0,num_v*10)
    for (k in 1:num_v) {
      for (l in neighb[[k]]) {
        i2[count] <- k
        j2[count] <- l
        z2[count] <-  -(num_neighb[[k]] +  num_neighb[[l]] -  sum(duplicated(c(neighb[[l]],neighb[[k]]))))
        count <- count + 1
      }
    }
    i2 <- i2[1:count-1]
    j2 <- j2[1:count-1]
    z2 <- z2[1:count-1]
    
    # Now add elements corresponding to the neighbours of the neighbours. Initialise, prob max 15 neighbours neighbours per node
    count <- 1
    i3 <- rep(0,num_v*15)
    j3 <- rep(0,num_v*15)
    z3 <- rep(0,num_v*15)
    neighb2 <- vector("list",num_v)
    for (k in 1:num_v) {
      # Construct neighbours of neighbours list, with redundancies (the value of the redundancy is then the element in Q)
      for (l in neighb[[k]]) {
        neighb2[[k]] <- c(neighb2[[k]],setdiff(neighb[[l]],c(neighb[[k]],k)))   # return number of elements in neighb[[l]] which are not in neighb[[k]] ( + node under consideration)
      }
      for (l in unique(neighb2[[k]]))  {
        i3[count] <- k
        j3[count] <- l
        z3[count] <- sum(neighb2[[k]] == l)
        count <- count + 1
      }
    }
    
    i3 <- i3[1:count-1]
    j3 <- j3[1:count-1]
    z3 <- z3[1:count-1]
    
    i <- c(i1,i2,i3)
    j <- c(j1,j2,j3)
    z <- c(z1,z2,z3)
    
  }
  
  
  
  z <- precinc*z
  Q <- sparseMatrix(i,j,x=z)
  return(Q)
}

# Create a covariance function from a set of locations
Covariance_from_points <- function(fn,pars,locs)  {
  n <- length(locs[[1]])
  V <- matrix(0,n,n)
  for(i in 1:n)
    for (j in 1:i)  {
      x1 = locs[[1]][i]
      y1 = locs[[2]][i]
      x2 = locs[[1]][j]
      y2 = locs[[2]][j]
      Dx <- x1-x2
      Dy <- y1-y2
      Dist <- sqrt(Dx^2 + Dy^2)
      if (fn == "sqexp") {
        V[j,i] <- V[i,j] <- pars[1]*exp(-Dist/pars[2])
      } else if (fn == "truncsqexp") {
        tau <- sqrt(2*pi/pars[2])
        if (tau*Dist > 2*pi) {
          V[j,i] <- V[i,j] <- 0
        } else {
          V[j,i] <- V[i,j] <- pars[1]/(3*pi)*((2*pi - tau*Dist)*(1 + (cos(tau*Dist))/2)  +
                                                3*sin(tau*Dist)/2)
        }
        
      }
      
    }
  return(V)
}

## Generate a random graph with num_v vertices with a mean of lambda neighbours each
Random_graph <- function(num_v,lambda) {
  neighb <- vector("list",num_v)
  for (i in c(1:num_v)) {
    num_neighb <- rpois(1,lambda = lambda)
    numattempts = 0
    while ((length(neighb[[i]]) <  num_neighb) && (numattempts < 10)) {
      elect_n <- floor(runif(1,min=i,max=num_v+1))
      numattempts <- numattempts + 1
      if (!(elect_n %in% neighb[[i]]) && i != elect_n) {
        neighb[[i]] <- c(neighb[[i]],elect_n)
        neighb[[elect_n]] <- c(neighb[[elect_n]],i)
      }
      
    }
  }
  return(neighb)
}


## Create a graph from a triangulation
Graph_from_tri <- function(p,tri) {
  
  neighb <- vector("list",dim(p)[1])
  for (i in 1:dim(tri)[1])  {
    this_tri = tri[i,]
    for (j in 1:3) {
      for (k in 1:3) {
        if (!(this_tri[k] %in% neighb[[this_tri[j]]]) && j != k)
          neighb[[this_tri[j]]] <- c(neighb[[this_tri[j]]],this_tri[k])
      }
    }
  }
  return(neighb)
}

## Create a graph from pixels
Graph_from_grid <- function(x,y,d) {
  n = length(x)
  X <- vector("list",n)
  for (i in 1:n) {
    neighb1 = intersect(which(x == x[i]),which(y == (y[i]+d)))
    neighb2 = intersect(which(x == x[i]),which(y == (y[i]-d)))
    neighb3 = intersect(which(x == (x[i]+d)),which(y == (y[i])))
    neighb4 = intersect(which(x == (x[i]-d)),which(y == (y[i])))
    X[[i]] <- c(neighb1,neighb2,neighb3,neighb4)
    
  }
  return(X)
}

## Explore a posterior function
explore_theta <- function(Covariance,log_theta_post,max_log_theta,dz=0.5,prune_grid = T) {
  # Standardize
  X <- eigen(Covariance,symmetric=T)
  ind <- sort(sort(diag(Covariance),decreasing=T,index.return=T)$ix,index.return=T)$ix
  X$values <- X$values[ind]
  X$vectors <- X$vectors[,ind]
  
  n <- dim(Covariance)[1]
  zexplore = matrix(rep(0,n))
  z_list <- vector("list",n)
  log_theta_post_max <- log_theta_post(max_par)
  if (n > 1) {
    VsqrtLambda <- X$vectors%*%sqrt(diag(X$values))
  } else {
    VsqrtLambda <- X$vectors%*%sqrt(X$values)
  }
  
  cat("Finding axial points...",sep="\n")
  flush.console()
  for (i in 1:n)    {
    while((log_theta_post_max - log_theta_post(c(max_par + VsqrtLambda%*%zexplore))) < 2.5) {
      z_list[[i]] <- c(z_list[[i]],zexplore[i])
      zexplore[i] <- zexplore[i] + dz
    }
    zexplore[i] = -dz
    while(((log_theta_post_max - log_theta_post(c(max_par + VsqrtLambda%*%zexplore))) < 2.5)
          &&(max_par[i] + (VsqrtLambda%*%zexplore)[i] > 0)) {
      z_list[[i]] <- c(z_list[[i]],zexplore[i])
      zexplore[i] <- zexplore[i] - dz
    }
    zexplore[i] <- 0
  }
  for (i in 1:n)    {
    z_list[[i]] <- sort(z_list[[i]])
  }
  
  cat("Forming grid...",sep="\n")
  flush.console()
  if (n == 1) {
    z = matrix(unlist(z_list))
  } else if (n == 2) {
    #Z <- meshgrid(z_list[[1]],z_list[[2]])
    z <- as.matrix(expand.grid(z_list[[1]],z_list[[2]]))
    # } else if (n == 3) {   #3d meshgrid NOT WORKING!!
    #    Z <- meshgrid(z_list[[1]],z_list[[2]],z_list[[3]])
    #    z <- matrix(c(c(Z$x),c(Z$y),c(Z$z)),length(Z$x),3)
  } else {
    nlist <- unlist(lapply(z_list,length))
    z <- matrix(0,prod(nlist),length(nlist))
    z[,length(nlist)] <- rep(z_list[[length(nlist)]],prod(nlist[1:(length(nlist)-1)]))
    for(i in seq(length(nlist),1,-1)) {
      if (i==1) {
        repmul = 1
      } else {
        repmul = prod(nlist[1:(i-1)])
      }
      temp <- rep(z_list[[i]],repmul)
      
      if ((i+1) > length(nlist)) {
        repmul = 1
      } else {
        repmul = prod(nlist[(i+1):length(nlist)])
      }
      
      whole_vec <- rep(temp,repmul)
      A <- reshape(as.matrix(whole_vec),length(temp),repmul)
      z[,i] = reshape(t(A),prod(nlist),1)
    }
  }
  
  theta_explore = t(max_par +  VsqrtLambda%*%t(z))
  d_theta <- VsqrtLambda%*%matrix(rep(dz,n))
  
  # Now we have a list of points on a grid and we can remove those which are not important
  cat("Computing log-posterior at all points and pruning...",sep="\n")
  flush.console()
  ind_keep <- NULL
  log_theta_post_vals <- rep(0,dim(theta_explore)[1])
  for (i in 1:length(theta_explore[,1])) {
    log_theta_post_vals[i] <- log_theta_post(theta_explore[i,])
    if ((log_theta_post_max - log_theta_post_vals[i]) < 2.5)
      ind_keep <- c(ind_keep,i)
    cat(paste(i/length(theta_explore[,1])*(100),"% complete"),sep="\n")
    flush.console()
  }
  
  
  if (prune_grid == T) {
    theta_explore <- theta_explore[ind_keep,]
    log_theta_post_vals <- log_theta_post_vals[ind_keep]
  }
  
  # Unnormalised posterior
  theta_post_unnorm <- exp(log_theta_post_vals - max(log_theta_post_vals))
  
  # Normalised posterior in theta scale
  vol <- abs(sum(theta_post_unnorm)*abs(prod(d_theta)))
  theta_post_norm <- theta_post_unnorm/vol
  
  # Plot marginals
  # Normalised posterior in z scals
  vol <- abs(sum(theta_post_unnorm)*abs(dz^n))
  theta_post_norm_z <- theta_post_unnorm/vol
  
  #for(i in 1:n) {
  # dev.new()
  # zaxis <- z[,i]
  # zmarg <- unlist(lapply(split(theta_post_norm_z,zaxis),function(x){sum(x*dz^(n-1)) }),use.names='F')
  # thetaaxis <- max_par[i] + VsqrtLambda[i,i]*unique(zaxis)
  # vol <- abs(sum(zmarg)*d_theta[i])
  # theta_marg_norm <- zmarg/vol
  # plot(thetaaxis,theta_marg_norm,type="o",xlab="theta_i",ylab="p(theta_i)")
  #}
  
  #theta_axis_all <- matrix(max_par,3,dim(z)[1]) + VsqrtLambda%*%t(z)
  thetaaxis <- vector("list",n)
  theta_marg_norm <- vector("list",n)
  for(i in 1:n) {
    #dev.new()
    zaxis <- z[,i]
    zmarg <- unlist(lapply(split(theta_post_norm_z,zaxis),function(x){sum(x*dz^(n-1)) }),use.names='F')
    
    zmat <- matrix(0,length(unique(zaxis)),n)
    zmat[,i] <- unique(zaxis)
    thetaaxis[[i]] <- (max_par + VsqrtLambda%*%t(zmat))[i,]
    vol <- abs(sum(zmarg)*d_theta[i])
    theta_marg_norm[[i]] <- zmarg/vol
    plot(thetaaxis[[i]],theta_marg_norm[[i]],type="o",xlab="theta_i",ylab="p(theta_i)")
  }
  
  
  return(list(post_unnorm = theta_post_unnorm,post_norm = theta_post_norm,
              theta_explore = theta_explore,d_theta = d_theta,
              thetaaxis = thetaaxis, theta_marg_norm = theta_marg_norm))
  
}

optimset <- function()
  return (list(maxit = 100L, tol=1e-7))

myoptim <- function(init,fn,gr,hessian,control=list(),aggressive=100,pars=NULL) {
  
  if (!is.null(pars)) {
    fn_wrap <- function(theta){ 
      return(fn(theta,pars))
    }
    gr_wrap <- function(theta){ 
      return(gr(theta,pars))
    }
  } else {
    gr_wrap = gr
    fn_wrap = fn
  }
  
  # Establish control parameter
  con <- optimset()
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  
  
  delta <- Inf
  theta <- init
  gradient <- gr_wrap(theta)
  thetanew <- init
  alpha = rep(0,length(theta))
  stepsize <- rep(0.01,length(theta))
  itercount <- 0
  while((abs(delta) > con$tol) && (stepsize > 1e-11) && (itercount < con$maxit)) {
    gradient <- gr_wrap(theta)
    for (i in 1:length(theta))   {
      alpha[i] <- gradient[i]
      # Find appropriate stepsize
      if (abs(alpha[i]) < stepsize[i])  { while (abs(alpha[i]) < stepsize[i]) alpha[i] <- alpha[i]*10
      } else if (abs(alpha[i])>(10*stepsize[i])) {while (abs(alpha[i])> (10* stepsize[i])) alpha[i] <- alpha[i]/10 }
      #thetanew[i] <- theta[i] - alpha[i]*abs(theta[i])
      thetanew[i] <- theta[i] - alpha[i]*stepsize[i]*aggressive
    }
    
    delta = as.vector(fn_wrap(thetanew) - fn_wrap(theta))
    if (delta < 0) {
      theta = thetanew
      stepsize <- rep(0.01,length(theta))
    } else {            # Overshoot detected
      stepsize <- stepsize/10
    }
    
    itercount <- itercount + 1
    cat(paste("Step ",itercount,": delta = ",delta,"theta = ",do.call(paste,as.list(theta))),sep="\n")
    flush.console()
  }
  
  if (abs(delta) < con$tol) cat("Reached local minimum",sep="\n")
  if (stepsize[1] < 1e-11) cat("Reached minimum stepsize (local minimum)",sep="\n")
  if (itercount == con$maxit) cat("Maximum iteration count reached",sep="\n")
  
  if (hessian == T) {
    H <- matrix(0,length(theta),length(theta))
    
    # get decent stepsizes in alpha (as above)
    stepsize <- 0.001
    for (i in 1:length(gradient)) {
      alpha[i] <- gradient[i]
      if (abs(alpha[i]) < stepsize)  { while (abs(alpha[i]) < stepsize) alpha[i] <- alpha[i]*10
      } else if (abs(alpha[i])>stepsize) {while (abs(alpha[i])> stepsize) alpha[i] <- alpha[i]/10 }
    }
    
    stepsizes <- abs(alpha*theta)
    
    
    for (i in 1:length(theta))
      for (j in 1:(i)) {
        maskvecj <- maskveci <- rep(0,length(theta))
        maskvecj[j] = stepsizes[j]
        maskveci[i] = stepsizes[i]
        gradient = 0.5*((gr_wrap(theta + .5*maskvecj)[i] - gr_wrap(theta - .5*maskvecj)[i])/(stepsizes[j]) +
                          (gr_wrap(theta + .5*maskveci)[j] - gr_wrap(theta - .5*maskveci)[j])/(stepsizes[i]) )
        H[j,i] <- H[i,j] <- gradient
        
      }
    
  } else { H = NA}
  return(list(par = theta,hessian = H))
}
marg_prec_from_kappa <- function(kappa_l,nu) {
  return(1/(gamma(nu)/(gamma(nu+1)*4*pi*kappa_l^(2*nu))))
}
Prec_from_SPDE <- function(M,K,tau,kappa,alpha=1)  {
  if (!(alpha %in% c(1,2,3,4))) {
    stop("alpha > 2 not implemented yet")
  }
  n <- nrow(M)
  if(class(tau) == "numeric") {
    tau <- sparseMatrix(i=1:n,j=1:n,x=tau)
  }
  if(class(kappa) == "numeric") {
    kappa <- sparseMatrix(i=1:n,j=1:n,x=kappa)
  }
  M_approx <- sparseMatrix(i=1:n,j=1:n,x=rowSums(M))
  M_approx_inv <- sparseMatrix(i=1:n,j=1:n,x=1/rowSums(M))
  M_kappa2 <- kappa%*%M%*%kappa
  G <- (M_kappa2 + K)
  if (alpha == 1) {
    Q <- tau%*%G%*%tau
  } else if (alpha ==2) {
    Q <- tau%*%G%*%M_approx_inv%*%G%*%tau
  } else if (alpha == 3) {
    Q <-  tau%*%G%*%M_approx_inv%*%G%*%M_approx_inv%*%G%*%tau
  } else if (alpha == 4) {
    Q <-  tau%*%G%*%M_approx_inv%*%G%*%M_approx_inv%*%G%*%M_approx_inv%*%G%*%tau
  }
  return(Q)
}

Build_AQA <- function(Qx,A,T) {
  # For now T >= 3
  
  n <- nrow(Qx)
  QA <- Qx %*% A
  AQA <- A %*% Qx %*% A + Qx
  AQ <- t(A) %*% Qx
  for ( i in 0:(T-3)) {
    if (i == 0) {
      Q <- cBind(-QA, AQA, -AQ , Zeromat(n,((T-3)-i)*n))
    } else if (i == (T-3)) {
      Q <- rBind(Q,cBind(Zeromat(n,n*i),-QA, AQA, -AQ))
    } else {
      Q <- rBind(Q,cBind(Zeromat(n,n*i),-QA, AQA, -AQ , Zeromat(n,((T-3)-i)*n)))
    }
  }
  Q <- #rBind(cBind(Qx, -AQ, Zeromat(n,(T-2)*n)),
        rBind(cBind((Qx - A %*% Qx %*% A) + A %*% Qx %*% A, -AQ, Zeromat(n,(T-2)*n)),
             Q,
             cBind(Zeromat(n,n*(T-2)),-QA, Qx))
  return(Q)
}
