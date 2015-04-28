#' @rdname dimensions
#' @aliases dimensions,measure-method
setMethod("dimensions",signature("measure"),function(obj){obj@dim})

#' @rdname dimensions
#' @aliases dimensions,manifold-method
setMethod("dimensions",signature("manifold"),function(obj){dimensions(obj@metric)})

#' @rdname dimensions
#' @aliases dimensions,domain-method
setMethod("dimensions",signature("domain"),function(obj){dimensions(obj@m)})


#' @rdname distance
#' @aliases distance,measure-method
setMethod("distance",signature("measure"),function(d,x1,x2){d@dist(x1,x2)})

#' @rdname distance
#' @aliases distance,manifold-method
setMethod("distance",signature("manifold"),function(d,x1,x2){distance(d@metric,x1,x2)})

# Distance function for polygon - taken from p_poly_dist (MATLAB)
mesh.dpoly <-  function(p,pv) {

  p1 <- p[,1]
  p2 <- p[,2]
  xv <- pv[,1]
  yv <- pv[,2]

  # If (xv,yv) is not closed, close it.
  xv <- c(xv);
  yv <- c(yv);
  Nv = length(xv);
  if ((xv[1] != xv[Nv]) || (yv[1] != yv[Nv])) {
    xv <- c(xv,xv[1])
    yv <- c(yv,yv[1])
    Nv <- Nv + 1;
  }

  # linear parameters of segments that connect the vertices
  A <- -diff(yv);
  B <-  diff(xv);
  C <- yv[-1]*xv[-Nv] - xv[-1]*yv[-Nv];

  dvec <- rep(0,length(p1))

  for (i in 1:length(p1)) {
    x <- p1[i]
    y <- p2[i]

    # find the projection of point (x,y) on each rib
    AB <- 1/(A^2 + B^2);
    vv <- (A*x+B*y+C);
    xp <- x - (A*AB)*vv;
    yp <- y - (B*AB)*vv;

    # find all cases where projected point is inside the segment
    idx_x <- (((xp>=xv[-Nv]) & (xp<=xv[-1])) | ((xp>=xv[-1]) & (xp<=xv[-Nv])));
    idx_y <- (((yp>=yv[-Nv]) & (yp<=yv[-1])) | ((yp>=yv[-1]) & (yp<=yv[-Nv])));
    idx <- idx_x & idx_y;

    # distance from point (x,y) to the vertices
    dv <- sqrt((xv[-Nv]-x)^2 + (yv[-Nv]-y)^2);

    if(!any(idx)) {# all projections are outside of polygon ribs
      d <- min(dv);
    } else {
      # distance from point (x,y) to the projection on ribs
      dp <- sqrt((xp[idx]-x)^2 + (yp[idx]-y)^2);
      d <- min(min(dv), min(dp));
    }


    if (in.poly(matrix(c(x,y),1,2),cbind(xv,yv))) {
      d = -d
    }
    dvec[i] <- d
  }
  return(matrix(dvec))

}
pol_from_xy <- function(x,y,order_vertices=T) {
  poly_list <- list()
  area <- NULL
  for (i in 1:nrow(x)) {
    this_edge <- t(rbind(x[i,],y[i,]))
    if(order_vertices) this_edge <- Order_vertices(this_edge)
    poly_list[[i]] <- as(this_edge,"gpc.poly")
  }
  return(poly_list=poly_list)
}
PolygonfromVoronoi <- function(Voronoi,p) {
  n = dim(p)[1]
  polygons <- vector("list",n)
  for (i in 1:n)  {
    #edges <- rbind(as.matrix(subset(Voronoi$dirsgs,ind1 ==i,select=c(x1,y1,bp1))),
    #         as.matrix(subset(Voronoi$dirsgs,ind2 ==i,select=c(x1,y1,bp1))),
    #         as.matrix(subset(Voronoi$dirsgs,ind1 ==i,select=c(x2,y2,bp2))),
    #         as.matrix(subset(Voronoi$dirsgs,ind2 ==i,select=c(x2,y2,bp2))))
    #X <- unique(edges)
    #if(sum(X[,3])>0) {
    #  X <- rbind(X,c(p[i,],0))
    #}
    #edges <- X[,(1:2)]

    X <- subset(Voronoi$dirsgs,(ind1 ==i | ind2 == i))
    X <- matrix(c(X$x1,X$x2,X$y1,X$y2,X$bp1,X$bp2),ncol=3)
    X <- unique(X)
    if(sum(X[,3])>0) {
      X <- rbind(X,c(p[i,],0))
    }
    edges <- X[,(1:2)]

    edges <- edges[chull(edges), ]
    polygons[[i]] <- as(edges,"gpc.poly")
  }
  return(polygons)

}
inpolygon_MATLAB <- function(x,y,xv,yv) {
  # vectorize the computation.
  # Translated from MATLAB
  inputSize = c(length(x),1);
  Nv <- length(xv)
  if (!any(is.nan(xv) | is.nan(yv))) {
    if ((xv[1] != xv[Nv]) || (yv[1] != yv[Nv])) {
      xv = c(xv,xv[1])
      yv = c(yv,yv[1])
      Nv = Nv + 1;
    }
  }

  mask = (x >= min(xv)) & (x <= max(xv)) & (y>=min(yv)) & (y<=max(yv));
  inbounds = which(mask);
  x = x[mask];
  y = y[mask];

  Np <- length(x)
  x <- t(matrix(x,Np,Nv))
  y <- t(matrix(y,Np,Nv))


  # Compute scale factors for eps that are based on the original vertex
  # locations. This ensures that the test points that lie on the boundary
  # will be evaluated using an appropriately scaled tolerance.
  # (m and mp1 will be reused for setting up adjacent vertices later on.)
  m = 1:(Nv-1)
  mp1 = 2:Nv
  avx = abs(0.5*(  xv[m] + xv[mp1]))
  avy = abs(0.5*(yv[m]+yv[mp1]))
  scaleFactor = pmax(avx[m], avy[m])
  scaleFactor = pmax(scaleFactor, avx[m]*avy[m])
  # Translate the vertices so that the test points are
  # at the origin.
  xv = matrix(xv,Nv,Np) - x
  yv = matrix(yv,Nv,Np) - y

  # Compute the quadrant number for the vertices relative
  # to the test points.
  posX = xv > 0
  posY = yv > 0
  negX = !posX
  negY = !posY
  quad = (negX & posY) + 2*(negX & negY) + 3*(posX & negY)

  # Ignore crossings between distinct edge loops that are separated by NaNs
  nanidx = which(is.nan(xv) | is.nan(yv));
  quad[nanidx] = NaN
  # Compute the sign() of the cross product and dot product
  # of adjacent vertices.
  theCrossProd = xv[m,]* yv[mp1,] - xv[mp1,]* yv[m,]
  signCrossProduct = sign(theCrossProd)


  # Adjust values that are within epsilon of the polygon boundary.
  # Making epsilon larger will treat points close to the boundary as
  # being "on" the boundary. A factor of 3 was found from experiment to be
  # a good margin to hedge against roundoff.
  scaledEps = scaleFactor*.Machine$double.eps*3
  idx = abs(theCrossProd) < scaledEps
  signCrossProduct[idx] = 0
  dotProduct = xv[m,]* xv[mp1,] + yv[m,]* yv[mp1,]

  # Compute the vertex quadrant changes for each test point.
  diffQuad = diff(quad)

  # Fix up the quadrant differences.  Replace 3 by -1 and -3 by 1.
  # Any quadrant difference with an absolute value of 2 should have
  # the same sign as the cross product.
  idx = (abs(diffQuad) == 3)
  diffQuad[idx] = -diffQuad[idx]/3
  idx = (abs(diffQuad) == 2)
  diffQuad[idx] = 2*signCrossProduct[idx]

  # Find the inside points.
  # Ignore crossings between distinct loops that are separated by NaNs
  nanidx = which(is.nan(diffQuad))
  diffQuad[nanidx] = 0
  inpoly = (apply(diffQuad,2,sum) != 0)

  # Find the points on the polygon.  If the cross product is 0 and
  # the dot product is nonpositive anywhere, then the corresponding
  # point must be on the contour.
  on = apply(((signCrossProduct == 0) & (dotProduct <= 0)),2,any)
  inpoly = inpoly | on

  mask[inbounds[!inpoly]] = 0
  return(c(reshape(matrix(mask),inputSize)))
}
inpolygon_MATLAB_big <- function(x,y,xv,yv) {
  # vectorize the computation.
  # Translated from MATLAB

  x_full <- x
  y_full <- y
  xv_orig <- xv
  yv_orig <- yv
  concat_dist <- NULL
  for (i in 1: ceil(length(x_full)/1000)) {
    this_batch <- ((i-1)*1000+1):min(c(i*1000,length(x_full)))
    x <- x_full[this_batch]
    y <- y_full[this_batch]
    xv <- xv_orig
    yv <- yv_orig

    inputSize = c(length(x),1);
    Nv <- length(xv)
    if (!any(is.nan(xv) | is.nan(yv))) {
      if ((xv[1] != xv[Nv]) || (yv[1] != yv[Nv])) {
        xv = c(xv,xv[1])
        yv = c(yv,yv[1])
        Nv = Nv + 1;
      }
    }

    mask = (x >= min(xv)) & (x <= max(xv)) & (y>=min(yv)) & (y<=max(yv));
    inbounds = which(mask);
    if( length(inbounds) == 0) {
      concat_dist <- rbind(concat_dist,matrix(0,length(this_batch),1))
    } else {
      x = x[mask];
      y = y[mask];

      Np <- length(x)
      x <- t(matrix(x,Np,Nv))
      y <- t(matrix(y,Np,Nv))


      # Compute scale factors for eps that are based on the original vertex
      # locations. This ensures that the test points that lie on the boundary
      # will be evaluated using an appropriately scaled tolerance.
      # (m and mp1 will be reused for setting up adjacent vertices later on.)
      m = 1:(Nv-1)
      mp1 = 2:Nv
      avx = abs(0.5*(  xv[m] + xv[mp1]))
      avy = abs(0.5*(yv[m]+yv[mp1]))
      scaleFactor = pmax(avx[m], avy[m])
      scaleFactor = pmax(scaleFactor, avx[m]*avy[m])
      # Translate the vertices so that the test points are
      # at the origin.
      xv = matrix(xv,Nv,Np) - x
      yv = matrix(yv,Nv,Np) - y

      # Compute the quadrant number for the vertices relative
      # to the test points.
      posX = xv > 0
      posY = yv > 0
      negX = !posX
      negY = !posY
      quad = (negX & posY) + 2*(negX & negY) + 3*(posX & negY)

      # Ignore crossings between distinct edge loops that are separated by NaNs
      nanidx = which(is.nan(xv) | is.nan(yv));
      quad[nanidx] = NaN
      # Compute the sign() of the cross product and dot product
      # of adjacent vertices.
      theCrossProd = xv[m,]* yv[mp1,] - xv[mp1,]* yv[m,]
      signCrossProduct = sign(theCrossProd)


      # Adjust values that are within epsilon of the polygon boundary.
      # Making epsilon larger will treat points close to the boundary as
      # being "on" the boundary. A factor of 3 was found from experiment to be
      # a good margin to hedge against roundoff.
      scaledEps = scaleFactor*.Machine$double.eps*3
      idx = abs(theCrossProd) < scaledEps
      signCrossProduct[idx] = 0
      dotProduct = xv[m,]* xv[mp1,] + yv[m,]* yv[mp1,]

      # Compute the vertex quadrant changes for each test point.
      diffQuad = diff(quad)

      # Fix up the quadrant differences.  Replace 3 by -1 and -3 by 1.
      # Any quadrant difference with an absolute value of 2 should have
      # the same sign as the cross product.
      idx = (abs(diffQuad) == 3)
      diffQuad[idx] = -diffQuad[idx]/3
      idx = (abs(diffQuad) == 2)
      diffQuad[idx] = 2*signCrossProduct[idx]

      # Find the inside points.
      # Ignore crossings between distinct loops that are separated by NaNs
      nanidx = which(is.nan(diffQuad))
      diffQuad[nanidx] = 0
      inpoly = (apply(diffQuad,2,sum) != 0)

      # Find the points on the polygon.  If the cross product is 0 and
      # the dot product is nonpositive anywhere, then the corresponding
      # point must be on the contour.
      on = apply(((signCrossProduct == 0) & (dotProduct <= 0)),2,any)
      inpoly = inpoly | on

      mask[inbounds[!inpoly]] = 0
      concat_dist <- rbind(concat_dist,reshape(matrix(mask),inputSize))
    }
  }
  return(c(concat_dist))
}

# Finds nearest neighbour of a continuous point from a set of discrete data points
interp_nn <- function(s,X) {
  s1 <- s[,1]
  s2 <- s[,2]
  x_all <- X[,1]
  y_all <- X[,2]
  z <- X[,3]

  # cater for nodes outside observation box
  outbox <- ((s1 < min(x_all)) | (s1 > max(x_all)) | (s2 < min(y_all)) |  (s2 > max(y_all)))
  if (length(s1)[1] == 1) {
    x <- z[which((abs(s1-x_all) == min(abs(s1-x_all)))&(abs(s2-y_all) == min(abs(s2-y_all))))[1]]  # Just in case it's in middle choose one corner
  } else {
    x <- apply(s,1,function(x) z[which((abs(x[1]-x_all) == min(abs(x[1]-x_all)))&
                                         (abs(x[2]-y_all) == min(abs(x[2]-y_all))))][1])
  }
  x[is.na(x)] <- 0
  return(x*!outbox)
}

# interpolates a continuous data set onto a grid
nn_grid_interp <- function(s,df,delta=10,miss_value=NA) {
  s_rnd <- data.frame(round(s/delta)*delta)
  s_rnd[,3] <- 1:nrow(s_rnd)
  names(s_rnd) = c("x","y","n")
  df_sub <- subset(df,(x >= min(s_rnd$x) &
                         x <= max(s_rnd$x) &
                         y >= min(s_rnd$y) &
                         y <= max(s_rnd$y)))
  df <- merge(df_sub,s_rnd,all.y=T)
  df <- arrange(df,n)
  if (any(is.na(df$z)))  df[is.na(df$z),]$z <- miss_value
  return(df$z)
}

# Regress gridded data on a triangulation
Regress_triang <- function(p,tri,z) {
  C_mean <- FindC(p,tri,z[,c(1,2)],method="C")
  EmptyC <-  which(colSums(C_mean) == 0)
  C_mean2 <- C_mean[,-EmptyC]
  n <- dim(p)[1]
  x_prior <- matrix(rep(0,n))
  x_prior[setdiff(1:n,EmptyC),] <- matrix(chol2inv(chol(t(C_mean2)%*%C_mean2))%*%t(C_mean2)%*%z[,3])
  return(x_prior)

  #C_mean <- FindC(p,tri,z[,c(1,2)],method="C")
  #EmptyC <-  which( apply(C_mean,2,sum) == 0)
  #z <- rbind(z,cbind(p[EmptyC,],0))
  #C_mean <- FindC(p,tri,z[,c(1,2)],method="C")
  #n <- dim(p)[1]
  #x_prior <- matrix(rep(0,n))
  #x_prior <- matrix(solve(t(C_mean)%*%C_mean)%*%t(C_mean)%*%z[,3])
  #return(x_prior)
}


# Check if there is at least one corner of a square inside a polygon
corner_in_poly <- function(cx,cy,lx,ly,poly) {
  return(which(!(inpolygon_MATLAB(cx+lx/2,cy+ly/2,poly[,1],poly[,2]) |  # Remove mascons out at sea
                   inpolygon_MATLAB(cx+lx/2,cy-ly/2,poly[,1],poly[,2]) |
                   inpolygon_MATLAB(cx-lx/2,cy+ly/2,poly[,1],poly[,2]) |
                   inpolygon_MATLAB(cx-lx/2,cy-ly/2,poly[,1],poly[,2]))))
}

# Check if polygons intersect
poly_in_poly <- function(pol_list,poly) {
  x <- rep(0,length(pol_list))
  for (i in 1:length(pol_list)) {
    x[i] <- area.poly(intersect(pol_list[[i]],poly)) - area.poly(pol_list[[i]])
  }
  return(which(x<0))
}

#Order vertices
Order_vertices <- function(p) {
  p2 <- p
  for (i in 1:(nrow(p)-1)) {



    dist_mat = as.matrix(dist(p2));
    if ((i > 3) && (dist_mat[1,i] < min(dist_mat[i,-c(1:i)])))  {
      break
    }
    else {
      ind = which(dist_mat[i,-(1:i)] == min(dist_mat[i,-(1:i)]))+i;
      temp =  p2[i+1,];
      p2[i+1,] = p2[ind[1],];
      p2[ind[1],] = temp;
    }
  }
  #p2 <- p2[1:(i-1),]
  return(p2)
}


poly_xy <- function(poly) {
  return(cbind(poly[[1]]@pts[[1]]$x,poly[[1]]@pts[[1]]$y))
}


# Piecewise interpolation of a function using tent functions.
# X is a data-frame with fields 'x' and 'y'
# ds denotes the spacing of the knots (1 by default)
# "first" denotes the location of the first knot (e.g. 2006)

Piecewise_interp <- function(X,ds=1,first=NULL,knots=NULL) {

  stopifnot(is(X,"data.frame"))
  stopifnot("x" %in% names(X) )
  stopifnot("y" %in% names(X) )
  stopifnot(is(ds,"numeric"))
  stopifnot(is.null(first) | is(first,"numeric"))
  stopifnot(is.null(knots) | is(knots,"numeric"))

  if (is.null(knots))
    if (is.null(first)) {             # Automatically place knots
      first_val <- ceiling(min(X$x)/ds)*ds - ds
      last_val <- floor(max(X$x)/ds)*ds + ds
      knots <- seq(first_val,last_val,by=ds)  # Form notes
    } else {    # Starting point pre-set
      first_val <- first
      knots <- first_val
      while( tail(knots,n=1) < max(X$x)) {
        knots <- c(knots,tail(knots,n=1)+ds)
      }
      last_val <- tail(knots,n=1)
      X <- subset(X,x > first_val)
    }
  # For each point, find basis which they belong to
  #X$section <- cut(X$x,knots)
  #count <- 1
  #X <- ddply(X,"section",function(df) {
  #  df$sectionnum = count
  #  count <<- count+1
  #  return(df)
  #}   )
  X$sectionnum <- cut(X$x,knots,labels=F)

  # If some sections are "unobserved", delete appropriately

  # For each section calculate values of phi
  X <- ddply(X,"section",function(df) {
    lo_x = knots[df$sectionnum]
    hi_x = knots[df$sectionnum+1]
    c1 = 1 - (-1/ds)*lo_x
    c2 = 1 - (1/ds)*hi_x

    df$phi1 <- (-1/ds)*df$x + c1
    df$phi1_ind <- df$sectionnum
    df$phi2 <- (1/ds)*df$x + c2
    df$phi2_ind <- df$sectionnum  + 1

    return(df)
  } )

  # Build the 'observation matrix'
  C <- sparseMatrix(i = c(1:nrow(X),1:nrow(X)), j = c(X$phi1_ind, X$phi2_ind), x = c(X$phi1, X$phi2))

  # Inference with EM algorithm for observation precision
  gamma_prec <- rep(1/var(X$y),10) #Initialise precision
  maxiter=100

  for (m in 2:maxiter) {
    # Estimate x
    ymat <- matrix(X$y)
    xML <- data.frame(x=knots,y=as.vector(chol2inv(chol(t(C) %*% C)) %*% t(C) %*% ymat))
    Qobs <- sparseMatrix(i=1:nrow(X),j=1:nrow(X),x=gamma_prec[m-1])
    Qtot <- t(C)%*%Qobs%*%C
    Varx <- solve(Qtot)

    # Update gamma_prec
    Gamma = as.vector(t(ymat) %*% ymat - 2*t(xML$y)%*% t(C) %*% ymat + sum(diag(t(C) %*% C %*% (Varx + xML$y%*%t(xML$y)))))
    gamma_prec[m] <- nrow(X)/Gamma
    if (abs(gamma_prec[m] - gamma_prec[m-1]) < 1e-5) break
  }
  cat(paste("EM converged after",m,"iterations",sep=" "),sep="\n")

  Vartrends <- NULL
  for (i in 1:(nrow(Varx)-1)) {
    Vartrends[i] <- (sum(diag(Varx)[i:(i+1)]) -2*(Varx[i,i+1]))/ds^2
  }
  trends = data.frame(year = head(knots,n=(length(knots)-1)),
                      trend = diff(xML$y)/ds,
                      std = sqrt(Vartrends)
  )
  return(list (x=xML, trends=trends))
  # returns a list with maximum likelihood estimate for the weights and the trend estimates
}

#' @title Find convex hull from a set of points
#' @author Andrew Zammit Mangion
#' @description Returns the convex hull from a set of points with fields \code{x} and \code{y}
#' @param points a data frame with fields \code{x} and \code{y}
#' @export
#' @return subset of \code{points} corresponding to the convex hull
#' @examples
#' N = 99
#' points <- data.frame(x = runif(n=N), y = runif(n=N),id=1:N)
#' hull <- find_hull(points)
#' plot(points$x,points$y)
#' lines(hull$x,hull$y,col="red")
find_hull <- function(points) {
  return(points[chull(points$x, points$y), ])
}

#' @title Find cpartition neighbourhood from underlying graph
#' @author Andrew Zammit Mangion
#' @description This function takes a set of points with a field \code{class} and the neighbourhood list of those points,
#' and computes the neighbourhood list of the super-graph composed by treating the \code{class} as a vertex in its own right. If even
#' one underlying vertex has a nerghbour in another partition, then that partition is the present partition's neighbour.
#' @param points a data frame with fields \code{x}, \code{y} and \code{class}
#' @param nei_list a list where each element is a vector containing the indices of the neighbours of the corresponding node
#' @export
#' @return a neighbourhood list of the superset
partition_nei_from_points <- function(points,nei_list) {
  partition_nei_list <- plyr::dlply(points,"class", function(df) {
    classes <- c()
    for (i in df$id) {
      points_to_examine <- nei_list[[i]]
      classes <- c(classes,points$class[points_to_examine])
    }
    classes <- classes[-which(classes == df$class[1])]
    return(unique(classes))
  })
  return(partition_nei_list)
}

#' @title Recursion for the Kernigham-Lin equations
#' @author C. Ladroue, A. Zammit Mangion
#' @description The  Kernigham-Lin equations split a graph into two equal parts depending on a weighting matrix. The weighting
#' matrix here is constrained to be the Eucledian distance between the points.
#' @param points a data frame with fields \code{x}, \code{y} and \code{id}, where \code{id} typically ranges from 1 to the number of points in the graph
#' @param nrep the number of times to split the graph
#' @return a list of vectors, where each vector contains the \code{id}s of the nodes in a particular group
#' @export
#' @examples
#' require(ggplot2)
#' N = 99
#' points <- data.frame(x = runif(n=N), y = runif(n=N),id=1:N)
#' points <- partition_section(points,nrep=4)
#' ggplot() + geom_point(data=points,aes(x,y,colour=class,size=4))
partition_section <- function(points,nrep) {
  partition <- .partition_section(points,nrep=nrep)
  points$class=0
  for( i in 1:(2^nrep)) {
    points$class[partition[[i]]] <- i
  }
  return(points)
}

#' @title Greedy graph colouring algorithm
#' @author Andrew Zammit Mangion
#' @description This algorithm colours the nodes in sequential order, where the order is given through a breadth first search algorithm started on a
#' pre-specified start-node.
#' @param adj.matrix the graphs adjacency matrix (usually sparse)
#' @param numcolours the maximum number of colours to try out with the \code{BFS} and \code{DSATUR} methods. The function gives an error if this is exceeded.
#' @param method takes values \code{"BFS"}, \code{"DSATUR"} or \code{"backtrack"} for the breadth-first-search, the maximum-degree-of-saturation and the backtracking respectively.
#' In the latter method, a four-colour configuration is attempted using brute force, where the order is obtained by first running a \code{DSATUR} run.
#' @param obs an \code{m x n} matrix identifiying which partitions each observations affect. If present, the algorithm
#' will not produce a colouring where an observation influences two or more partitions of the same colour
#' (as such this increases the chromatic number of the problem). This option is only available with \code{BFS} and
#' \code{DSATUR}
#' @param startnode the starting node for the BFS algorithm
#' @return a data frame with the associated colour for each vertex (also \code{class})
#' @references \url{http://community.topcoder.com/longcontest/?module=Static&d1=match_editorials&d2=intel_mtcs_10}
#' @export
colour_graph <- function(adj.matrix,numcolours,method="BFS",startnode=1,obs=NULL) {
  partition_colour <- data.frame(class = 1:nrow(adj.matrix),colour=0)
  G <- graph.adjacency(adj.matrix)

  if(method == "BFS") {
    G.bfs <- graph.bfs(G,startnode)
    for(i in G.bfs$order) {
      # check neighbour colours
      # use minimum of intersection of 1:numcolours and neighbour colours as this colour
      nei_colours <- partition_colour$colour[partition_nei_list[[i]]]
      if (suppressWarnings(any(nei_colours))) {
        if(!(is.null(obs))) {
          offending_obs <- which(obs[,i] == 1)
          affected_partitions <- Reduce("union",apply(D[offending_obs,],1,function(x) which(x==1)))
          prohibited_colours <- partition_colour$colour[affected_partitions]
          partition_colour$colour[i] <- min(setdiff(c(1:numcolours),union(nei_colours,prohibited_colours)))
        } else {
          partition_colour$colour[i] <- min(setdiff(c(1:numcolours),nei_colours))
        }
      } else {
        partition_colour$colour[i] <- 1
      }
    }
  }


  if(method %in% c("DSATUR","backtrack")) {
    degrees <- rowSums(adj.matrix)
    i <- which.max(degrees)
    col.matrix <- adj.matrix
    done_list <- NULL
    while(any(partition_colour$colour == 0)) {
      done_list <- c(done_list,i)
      nei_colours <- partition_colour$colour[partition_nei_list[[i]]]
      ## Color node
      if (suppressWarnings(any(nei_colours))) {
        if(!(is.null(obs))) {
          offending_obs <- which(D[,i] == 1)
          affected_partitions <- Reduce("union",apply(D[offending_obs,],1,function(x) which(x==1)))
          prohibited_colours <- partition_colour$colour[affected_partitions]
          partition_colour$colour[i] <- min(setdiff(c(1:numcolours),union(nei_colours,prohibited_colours)))
        } else {
          partition_colour$colour[i] <- min(setdiff(c(1:numcolours),nei_colours))
        }
      } else {
        partition_colour$colour[i] <- 1
      }

      if(length(done_list) < nrow(adj.matrix)) {

        col.matrix[which(col.matrix[,i]>0),i] <-  partition_colour$colour[i]+1
        ## Find next node
        sat.level <- apply(col.matrix,1,function(x) length(unique(x)))
        maxsat <- max(sat.level)
        found_node <- 0
        while(!found_node) {
          vv <- setdiff(which(sat.level == maxsat),done_list)
          if(length(vv) == 0) {
            maxsat <- maxsat - 1
          } else {
            i <- vv[which.max(degrees[vv])]
            found_node <- 1
          }
        }
      }
    }
  }

  if(method=="backtrack") {
    numcolours=4
    partition_colour <- data.frame(class = 1:nrow(adj.matrix),colour=0)
    order <- done_list
    done <- 0
    i <- order[1]
    ilist <- NULL
    tried_colours <- list(); length(tried_colours) <- nrow(adj.matrix)
    program_state <- list()
    global_count <- 0
    while (!done) {
      backtrack = 0
      nei_colours <- partition_colour$colour[partition_nei_list[[i]]]
      if (suppressWarnings(any(nei_colours))) {
        colours_to_try <- setdiff(1:numcolours,tried_colours[[i]])
        if (length(setdiff(colours_to_try,nei_colours)) ==0) {
          ## BACKTRACK
          X <- match(partition_nei_list[[i]],ilist) # return first positions of neighbours in ilist
          go.back.to <- max(X,na.rm=T)
          i <- ilist[go.back.to]
          #restore state
          tried_colours <-  program_state[[go.back.to]]$tried_colours
          partition_colour <-  program_state[[go.back.to]]$partition_colour
          global_count <- go.back.to - 1
          backtrack = 1
          ilist <- ilist[1:go.back.to-1] ##### TEST THIS LINE
        } else {
          colour <- min(setdiff(colours_to_try,nei_colours))
          if(global_count + 1 == nrow(adj.matrix)) done = 1
        }

      } else {
        colour <- min(setdiff(c(1 :numcolours),tried_colours[[i]]))
      }

      if(!backtrack) {
        global_count <- global_count + 1
        ilist <- c(ilist,i)
        partition_colour$colour[i] <- colour
        tried_colours[[i]] <- c(tried_colours[[i]],colour)
        i <- order[which(order == i) + 1]
        program_state[[global_count]] <- list(tried_colours = tried_colours, partition_colour = partition_colour)
      }
      print(i)
    }
  }


  return(partition_colour)
}


.partition_section <- function(points,nrep) {
  made_odd = 0
  if(is.odd(nrow(points))) {  # we need to make it even
    dummy_id <- tail(sort(points$id),1) + 1
    points <- rbind.fill(points,data.frame(x = mean(points$x), y = mean(points$y),id=dummy_id))
    made_odd = 1
  }

  weightMatrix <- as.matrix(dist(subset(points,select=c("x","y"))))
  partition <- approximateBisection(-weightMatrix)

  if(made_odd) {
    in_part <- any(points[partition[[1]],]$id  == dummy_id)*1 + any(points[partition[[2]],]$id  == dummy_id)*2
    ind_to_remove <- which(points[partition[[in_part]],]$id  == dummy_id)
    partition[[in_part]] <- partition[[in_part]][-ind_to_remove]
  }

  if (nrep == 1) {
    return(list(points$id[partition[[1]]],points$id[partition[[2]]]))
  } else {
    return(c(.partition_section(points[partition[[1]],],nrep = nrep - 1),
             .partition_section(points[partition[[2]],],nrep = nrep - 1)))
  }
}

# Approximate bisection
# returns a bisection of the graph that minimizes the cost using Kernighan/Lin Algorithm
# http://www.eecs.berkeley.edu/~demmel/cs267/lecture18/lecture18.html#link_4.2
# partition<-approximateBisection(weightMatrix)
# weightMatrix is symmetric matrix of size 2Nx2N made of non-negative values.
# partition is a list of two vectors of N indices.
# C.Ladroue
approximateBisection<-function(weightMatrix,mode="matrix",minimumGain=1e-5){
  #   minimumGain<-1e-5 # minimum value for gain, setting it to 0 might lead to infinite loop due to numerical inaccuracy

  N<-dim(weightMatrix)[1] # number of elements
  m<-N/2

  # start off with a random partition
  A<-sample(1:N,N/2,replace=FALSE)
  B<-(1:N)[-A]

  maxGain<-Inf
  while(maxGain>minimumGain){
    DA<-rowSums(weightMatrix[A,B])-rowSums(weightMatrix[A,A])+diag(weightMatrix[A,A])
    DB<-rowSums(weightMatrix[B,A])-rowSums(weightMatrix[B,B])+diag(weightMatrix[B,B])
    unmarkedA<-1:m
    unmarkedB<-1:m
    markedA<-rep(0,m)
    markedB<-rep(0,m)
    gains<-rep(0,m)
    for(k in 1:m){
      # find best pair from remainder
      # with 2 loops, slow but easy on memory
      if(mode=='2loops'){
        bestGain<--Inf
        besti<-0
        bestj<-0
        for(i in unmarkedA)
          for(j in unmarkedB){
            gain<-DA[i]+DB[j]-2*weightMatrix[A[i],B[j]]
            if(gain>bestGain) {bestGain<-gain; besti<-i;bestj<-j}
          }
        #           mark the best pair
        unmarkedA<-unmarkedA[-which(unmarkedA==besti)]
        unmarkedB<-unmarkedB[-which(unmarkedB==bestj)]
        markedA[k]<-besti
        markedB[k]<-bestj
      }
      # with one matrix, much faster but builds a matrix as large as weightMatrix
      if(mode=='matrix'){
        dimension<-m+1-k
        fasterGain<-matrix(DA[unmarkedA],nrow=dimension,ncol=dimension,byrow=FALSE)+
          matrix(DB[unmarkedB],nrow=dimension,ncol=dimension,byrow=TRUE)-
          2*weightMatrix[A[unmarkedA],B[unmarkedB]]
        # mark the best pair
        best<-arrayInd(which.max(fasterGain),.dim=c(dimension,dimension))
        besti<-unmarkedA[best[1]]
        bestj<-unmarkedB[best[2]]
        bestGain<-fasterGain[best]
        markedA[k]<-unmarkedA[best[1]]
        markedB[k]<-unmarkedB[best[2]]
        unmarkedA<-unmarkedA[-best[1]]
        unmarkedB<-unmarkedB[-best[2]]
      }
      # record gain
      gains[k]<-bestGain

      # update D for unmarked indices
      DA[unmarkedA]<-DA[unmarkedA]+2*weightMatrix[A[unmarkedA],A[besti]]-2*weightMatrix[A[unmarkedA],B[bestj]]
      DB[unmarkedB]<-DB[unmarkedB]+2*weightMatrix[B[unmarkedB],B[bestj]]-2*weightMatrix[B[unmarkedB],A[besti]]
    }
    gains<-cumsum(gains)
    bestPartition<-which.max(gains)
    maxGain<-gains[bestPartition]
    if(maxGain>minimumGain){
      # swap best pairs
      A1<-c(A[-markedA[1:bestPartition]],B[markedB[1:bestPartition]])
      B1<-c(B[-markedB[1:bestPartition]],A[markedA[1:bestPartition]])
      A<-A1
      B<-B1
    }
  }
  list(A,B)
}


#' @title Extract border
#' @author Andrew Zammit Mangion
#' @description The border is extracted from shapes supplied in x-y format
#' @param df a data frame with fields \code{x} and \code{y}.
#' @param delta the grid spacing.
#' @details The algorithm works by finding how many neighbours a grid point has using \code{fields::fields.rdist.near}. A grid point with less than 4 neighbours is assumed to be a border point. For this algorithm to work, it is important that the shape is `full' and not contain any missing values.
#' @return a subset of \code{df} containing the border points.
#' @export
#' @examples
#' grid <- expand.grid(x = c(1:10), y = c(1:10))
#' names(grid) = c("x","y")
#' border_points <- border_from_grid_points(df=grid,delta=1)
border_from_grid_points <- function(df,delta) {
  D <- fields.rdist.near(df[c("x","y")],df[c("x","y")],delta=delta)  # Find distance matrix
  D2 <- sparseMatrix(i=D$ind[,1],j=D$ind[,2],x=1)                        # See how many neighbours
  nn <- rowSums(D2,1,sum)
  border_points <- df[which(nn < 5),]
}
