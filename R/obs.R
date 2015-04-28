#' @title Pre-processes observations
#' @description This function simply takes a data frame and performs standard pre-processing functions such as removing obvious outliers and averaging data
#' in computationally efficient way over a regular grid
#'
#' @param Obs.obj The observation object
#' @param std_lim Values with error higher than \code{std_lim} are deleted
#' @param abs_lim Values which are bigger than \code{abs_lim} are deleted
#' @param avr_method If \code{mean} then the mean value and mean error on a sub-grid are used to sub-sample the observations. If \code{median} the median value and MAD of the z-values are used
#' for subsampling. If \code{mean_no_std} or \code{median_no_std}, the error is not computed ignored.
#' @param box_size The grid width over which the observations are averaged
#' @param min_pts If a grid box contains less than \code{min_pts} it is marked as empty
#' @param ... Other parameters which are ignored
#' @export
preprocess_obs <- function(Obs.obj,std_lim=NA,abs_lim=NA,avr_method=NA,box_size=10,min_pts=4,...) {
  
  if(!(is.na(std_lim))) {
    stopifnot(is.numeric(std_lim))
    Obs.obj@df <- subset(Obs.obj@df,std < std_lim)
  }
  
  if (!is.na(abs_lim)) {
    stopifnot(is.numeric(abs_lim))
    Obs.obj@df <- subset(Obs.obj@df,abs(z) < abs_lim)
  }
  
  if(any(Obs.obj@df$std == 0)) {
    warning("Data points with zero error detected. These are automatically deleted")
    Obs.obj@df <- subset(Obs.obj@df,std > 0)
  }
  
  
  if(!is.na(avr_method)) {
    stopifnot(avr_method %in% c("mean","median","mean_no_std","median_no_std"))
    
    breaksx <- seq(min(Obs.obj@df$x)-1,max(Obs.obj@df$x)+1,by=box_size)
    breaksy <- seq(min(Obs.obj@df$y)-1,max(Obs.obj@df$y)+1,by=box_size)
    Obs.obj@df$box_x <- cut(Obs.obj@df$x,breaksx,labels=F)
    Obs.obj@df$box_y <- cut(Obs.obj@df$y,breaksy,labels=F) 
    
    
    averaging_fun <- function(z) {
      x <- z[1]  
      if(length(z) > min_pts) {
        if(avr_method %in% c("median","median_no_std")) {
          x <- median(z)
        } else if(avr_method  %in% c("mean","mean_no_std") ) {
          x <- mean(z)
        }
      }
      return(x)
    }
    
    std_fun <- function(z,std) {
      x <- std[1]  
      if(length(z) > 4) {
        if(avr_method=="median") {
          x <- mad(z)
        } else if(avr_method=="mean") {
          x <- sd(z)
        }
      }
      return(x)
    }
    
    boxes <- group_by(Obs.obj@df, box_x,box_y,t)
    if(avr_method %in% c("median","mean")) { 
      Obs.obj@df <- data.frame(summarise(boxes,x=round(mean(x)),y=round(mean(y)),z2=averaging_fun(z),std=std_fun(z,std)))
    } else if(avr_method %in% c("median_no_std","mean_no_std")) {
      Obs.obj@df <- data.frame(summarise(boxes,x=round(mean(x)),y=round(mean(y)),z2=averaging_fun(z)))
    }
    Obs.obj@df$z <- Obs.obj@df$z2
    Obs.obj@df$z2 <- NULL
    
    Obs.obj@df <- arrange(Obs.obj@df,t)
    Obs.obj@n <- 1:nrow(Obs.obj@df)
  }  
  return(Obs.obj)
}

#' @rdname setalpha
#' @aliases setalpha,Obs-method
setMethod("setalpha",signature(.Object="Obs"),
          function(.Object,alpha,av_dist) {
            .Object@args$alpha0 <- alpha
            .Object@args$av_dist <- av_dist
            .Object@args$P <-  Find_Smooth_mat(subset(.Object@df,t==.Object@df$t[1]),alpha,av_dist)
            return(.Object)          })


#' @rdname split_validation
#' @aliases split_validation,Obs-method
setMethod("split_validation",signature = "Obs",function(.Object,samples,common=0,...) {
  y_tot <- getDf(.Object)
  y_tot$n <- 1:nrow(y_tot)
  y_tot_sub <- subset(y_tot,...)
  
  if(common) {
    samples_per_t <- round(samples/length(unique(y_tot$t)))
    data_per_t <- nrow(subset(y_tot,t==0))
    XY <- unique(y_tot[c("x","y")])
    XY$n <- 1:nrow(XY)
    xylocs <- sample(XY$n,samples_per_t,replace=F)
    obs_test_indices <- c(daply(y_tot,"t",function(df) {
      return(df$n[xylocs])}))
  } else {
    obs_test_indices <- sample(y_tot_sub$n,samples,replace=F)
  }
  
  y_test <- y_tot[obs_test_indices,]
  y_test <- arrange(y_test,n)
  y_tot <- y_tot[-obs_test_indices,]
  y_tot <- arrange(y_tot,n)
  
  if(is(.Object,"Obs_poly")) {
    O_test <- .Object
    O_test@pol <- O_test@pol[obs_test_indices]
    O_test@df2 <- subset(O_test@df2,id %in% obs_test_indices)
    O_test@df <- subset(O_test@df,id %in% obs_test_indices)
    O_test@n <- length(obs_test_indices)
    O_test@args$P <-  Imat(O_test@n)
    
    O_pruned <- .Object
    O_pruned@pol <- O_pruned@pol[-obs_test_indices]
    O_pruned@df2 <- subset(O_pruned@df2,!(id %in% obs_test_indices))
    O_pruned@df <- subset(O_pruned@df,!(id %in% obs_test_indices))
    O_pruned@n <- length(nrow(y_tot) - obs_test_indices)
    O_pruned@args$P <-  Imat(O_pruned@n)
    
  } else {
    O_test <- new("Obs",df=y_test)
    O_pruned <- new("Obs",df=y_tot)
  }
  
  return(list(O_pruned = O_pruned, O_val = O_test))
  
})



#' @title Load data
#' @description Loads data simultaneously from multiple files. For data files, it is assumed that the columns have a header. The function understands filenames with extensions 'shp', 'txt', 'dat', and 'tif'.
#' @param thinning determines a thinning rate for the shapefiles. This is useful when shapefiles, particularly coastlines, are very high in resolution.
#' @param tif_grid the native resolution of the tif file, if present, in m
#' @param convert_m_to_km if \code{TRUE}, shapefiles and tif file units are converted to km
#' @param ... paths to data sets which will be loaded. These need to be assigned variable names
#' @return list of data frames
#' @export
#' @examples
#' df <- data.frame(x=c(1,2,3),y=c(2,2,2))
#' write.table(df,file="~/temp.dat")
#' data <- Load_data("~/temp.dat")
#' file.remove("~/temp.dat")
Load_data <- function(thinning = 1,     # how much we thin grounding and shapefiles by
                      tif_grid = 750,     # grid spacing of tif image in metres
                      convert_m_to_km = T, # works depending on extension (not if data is just a table)
                      ...) {
  
  args <- list(...)
  shapefiles <- list()
  
  for (i in seq_along(args)) {
    if (!file.exists(args[[i]])) stop(paste("Cannot find file ",args[[i]],sep=""))
    if (basename_ext(args[[i]]) == "shp") {
      cat("Loading shapefile...",sep="\n")
      shpfile <- readShapeSpatial(args[[i]])
      shpfile_table <- fortify(shpfile)
      if (convert_m_to_km) {
        shpfile_table$long <- shpfile_table$long/1000
        shpfile_table$lat <- shpfile_table$lat/1000
      }
      shpfile_table <- shpfile_table[seq(1,nrow(shpfile_table),thinning),]
      shp_names <- names(shpfile_table)
      shp_names[which(shp_names == "long")] <- "x"
      shp_names[which(shp_names == "lat")] <- "y"
      names(shpfile_table) <- shp_names
    } else if(basename_ext(args[[i]]) %in% c("txt","dat")) {
      cat("Loading table from file...",sep="\n")
      shpfile_table <- read.table(args[[i]],header=T)
    } else if(basename_ext(args[[i]]) %in% "tif") {
      cat("Loading tif file...",sep="\n")
      cat(paste("...Setting tif grid to ",tif_grid," m",sep=""),sep="\n")
      X <- raster(args[[i]])
      xmin <- X@extent@xmin
      xmax <- X@extent@xmax
      ymin <- X@extent@ymin
      ymax <- X@extent@ymax
      
      if (convert_m_to_km) {
        xmin <- xmin/1000
        xmax <- xmax/1000
        ymin <- ymin/1000
        ymax <- ymax/1000
        tif_grid <- tif_grid/1000
      }
      
      xgrid <- seq(xmin,by=tif_grid,length=X@ncols)
      ygrid <- seq(ymax,by=-tif_grid,length=X@nrows)
      shpfile_table <- expand.grid(xgrid,ygrid)
      names(shpfile_table) <- c("x","y")
      shpfile_table$z <- X[] 
    }
    shapefiles[[i]] <- shpfile_table
    names(shapefiles)[i] <- names(args)[i]
  }
  return(shapefiles)
  
}

#' Extract roughness from topography
#'
#' Supplied with a data frame with fields x, y and z, this function bins the data into boxes of length \code{ds} and returns a data frame where \code{z} is now the estimated roughness, given by the lof of the standard deviation of the topography in the boxes.
#' @param df data frame with columns \code{x, y} and \code{z}
#' @param ds the grid size required
#' @export
#' @examples
#' df <- data.frame(x=runif(1000),y=runif(1000),z=runif(1000))
#' df_roughness <- roughness_from_topo(df,ds=0.1)
roughness_from_topo <- function(df,ds=20) {
  
  
  names(df)[ncol(df)] <- "topo_el"
  xmin <- min(df$x)
  xmax <- max(df$x)
  ymin <- min(df$y)
  ymax <- max(df$y)
  
  first_val <- ceiling(xmin/ds)*ds - ds
  last_val <- floor(xmax/ds)*ds + ds
  x_super_grid <- seq(first_val,last_val,by=ds)
  
  first_val <- ceiling(ymin/ds)*ds - ds
  last_val <- floor(ymax/ds)*ds + ds
  y_super_grid <- seq(first_val,last_val,by=ds)
  
  df$box_x <- cut(df$x,x_super_grid,labels=F)
  df$box_y <- cut(df$y,y_super_grid,labels=F)
  
  
  f1 <- function(x) {
    if(length(x) > 4){
      x <- ceiling(x[1]/ds)*ds - ds 
    } else { x = x[1] }
  }
  
  f2 <- function(x) {
    if(length(x) > 4){
      x <- log(sd(x) + 1)
    } else { x = NA }
  }
  
  boxes <- group_by(df, box_x,box_y)
  df <- data.frame(summarise(boxes,x=f1(x),y=f1(y),z=f2(topo_el)))
  df <- df[order(df$x,df$y),]
  
  
  return(df)
}



Find_Smooth_mat <- function(df,alpha0,av_dist=450){
  m <- nrow(df)
  if(alpha0 == 0) {
    return(Imat(m)) 
  } else {
    D <- fields::fields.rdist.near(cbind(df$x,df$y),cbind(df$x,df$y),delta=av_dist)
    X2 <- sparseMatrix(i=D$ind[,1],j=D$ind[,2],x=D$ra) 
    P <- sparseMatrix(i=D$ind[,1],j=D$ind[,2],x=alpha0) 
    diag(P) <- apply(P,1,function(x) 1 - sum(x) + alpha0)
    P <- as(P,"dgTMatrix")
    dense_obs <- which(diag(P) < 0.1)
    for (i in dense_obs) {
      indices <- which((P@i + 1) ==i)
      P@x[indices] <- 1/length(indices)
    }
  }
  return(P)}


.expand_poly <- function(df) {
  df <- df[,!(names(df) %in% c("x","y"))]
  x_ind <- grep("x[[:digit:]]",names(df))
  y_ind <- grep("y[[:digit:]]",names(df))
  df3 <- adply(df,1,function(dd) {xy <- data.frame(x = as.numeric(dd[,x_ind]),y=as.numeric(dd[,y_ind]))
                                  dd <- dd[,-c(x_ind,y_ind)]
                                  row.names(dd) <- row.names(xy) <- NULL
                                  return(cbind(dd,xy)) })
  x_ind <- grep("x[[:digit:]]",names(df3))
  y_ind <- grep("y[[:digit:]]",names(df3))
  return(df3[,-c(x_ind,y_ind)])
}
