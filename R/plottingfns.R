# FRK: An R Software package for spatial and spatio-temporal prediction
# with large datasets.
# Copyright (c) 2017 University of Wollongong
# Author: Andrew Zammit-Mangion, azm (at) uow.edu.au
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#' @title Draw a map of the world with country boundaries.
#' @description Layers a \code{ggplot2} map of the world over the current \code{ggplot2} object.
#' @param g initial ggplot object
#' @param inc_border flag indicating whether a map border should be drawn or not; see details.
#' @details This function uses \code{ggplot2::map_data} in order to create a world map. Since, by default, this creates lines crossing the world at the (-180,180) longitude boundary, function \code{.homogenise_maps} is used to split the polygons at this boundary into two. If \code{inc_border} is TRUE, then a border is drawn around the lon-lat space; this option is most useful for projections that do not yield rectangular plots (e.g., the sinusoidal global projection).
#' @seealso the help file for the dataset \code{\link{worldmap}}
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' draw_world(g = ggplot())}
draw_world <- function(g = ggplot() + theme_bw() + xlab("") + ylab(""),inc_border = TRUE) {

    ## Basic checks
    if(!(is(g, "ggplot"))) stop("g has to be of class ggplot")
    if(!(is.logical(inc_border))) stop("inc_border needs to be TRUE or FALSE")

    ## Suppress bindings warning
    long <- lat <- group <- NULL

    ## Load the world map data from the FRK package
    data(worldmap, envir=environment(), package = "FRK")

    ## Homogenise (see details) to avoid lines crossing the map
    worldmap <- .homogenise_maps(worldmap)

    ## If user wants to draw border
    if(inc_border) {

        ## Create a border data frame at lon/lat boundaries
        border <- data.frame(long=c(-179.99,-179.99,179.99,179.99),
                             lat=c(-89.99,89.99,89.99,-89.99),
                             group=1e5,       # create a high group number to avoid duplication
                             region="border") # create a new name for border not already used
        worldmap <- plyr::rbind.fill(worldmap,border) # just append it to world map
    }

    ## Now return a gg object with the map overlayed
    g + geom_path(data = worldmap, aes(x=long, y=lat, group=group), colour="black",size=0.1)
}

#' @rdname show_basis
#' @aliases show_basis,Basis-method
setMethod("show_basis",signature(basis = "Basis"),  # GRBF basis with mean offset as last weight
          function(basis,g=ggplot() + theme_bw() + xlab("") + ylab("")) {

              ## Currently only spherical basis functions are plotted. In principle, the manifold might
              ## be changed to reflect anisotropy/heterogeneity and the below plotting functions
              ## Suppress bindings warning
              y <- res <- x <- lon <- lat <- NULL

              ## If we are on the real line
              if(is(manifold(basis),"real_line")) {
                  s1min <- min(basis@df$loc1)  - max(basis@df$scale)*3  # suitable minimum of s
                  s1max <- max(basis@df$loc1)  + max(basis@df$scale)*3  # suitable maximum of s
                  s <- matrix(seq(s1min,s1max,length=1000))             # create s-axis
                  for (i in 1:basis@n) {                                # for each basis function
                      S <- basis@fn[[i]](s)                             # evaluate fuction over s-axis
                      df <- data.frame(s=as.numeric(s),                 # create data frame with
                                       y = as.numeric(S),               # basis function
                                       res=basis@df$res[i])

                      ## Draw gg object
                      g <- g + geom_line(data=df,aes(x=s,y=y,col=as.factor(res))) +
                          labs(colour="res")
                  }

                  ## If we are on the plane
              } else  if(is(manifold(basis),"plane")) {

                  ## can be amended eventually to reflect anisotropy etc.
                  message("Note: show_basis assumes spherical distance functions when plotting")

                  l <- lapply(1:basis@n,function(i) {   # for each basis function

                      ## Create a data frame containin the x,y coordinates of a circle
                      ## around the basis function centroid and the function's resolution
                      data.frame(.circleFun(center=as.numeric(basis@df[i,1:2]),
                                            diameter = basis@df$scale[i]),
                                 res=basis@df$res[i],
                                 id = i)})
                  df <- bind_rows(l)              # quick rbind of l
                  df$res <- as.factor(df$res)     # convert to factor

                  ## Draw circles with different linetypes for the different resolutions
                  g <- g + geom_path(data=df,
                                     aes(x=x,y=y,group=id,linetype=res))

              } else  if(is(manifold(basis),"sphere")) {
                  ## If we're on the sphere we just show circles proportional in size to the resolution as
                  ## it makes for a neater figure
                  df <-data.frame(basis)                    # extract data frame
                  df <- df[rev(rownames(df)),]              # reverse order of data frame
                  names(df)[1:2] <- c("lon","lat")          # ensure the first two columns are labelled correctly

                  ## Draw the circles in lon and lat
                  g <- g + geom_point(data=df,aes(x=lon,y=lat,size=res),shape=1) +
                      scale_size_continuous(trans="reverse",breaks =1:10)


                  ## If we're on the space-time plane do as above but draw the bases at each time point
                  ## Note: This is never used as we always have Tensor Basis in practice (see below)
              } else  if(is(manifold(basis),"STplane")) {
                  df <-basis@df
                  df <- df[rev(rownames(df)),]
                  names(df)[1:2] <- c("x","y")                # ensure the first two columns are labelled correctly
                  g <- g + geom_point(data=df,aes(x=x,y=y,size=res),shape=1) +
                      scale_size_continuous(trans="reverse",breaks =1:10) +
                      facet_wrap(~loc3)

                  ## If we're on the space-time plane draw the bases at each time point
                  ## Note: This is never used as we always have Tensor Basis in practice (see below)
              } else  if(is(manifold(basis),"STsphere")) {
                  df <-basis@df
                  df <- df[rev(rownames(df)),]
                  names(df)[1:2] <- c("lon","lat")            # ensure the first two columns are labelled correctly
                  g <- g + geom_point(data=df,aes(x=lon,y=lat,size=res),shape=1) +
                      scale_size_continuous(trans="reverse",breaks =1:10) +
                      facet_wrap(~loc3)

              }
              return(g + theme_bw())

          })

#' @rdname show_basis
#' @aliases show_basis,TensorP_Basis-method
setMethod("show_basis",signature(basis = "TensorP_Basis"),
          function(basis,g=ggplot()) {
              ## For Tensor Basis just plot first the spatial and then the temporal
              (show_basis(basis@Basis1) + ggtitle("Basis1")) %>% print()
              (show_basis(basis@Basis2) + ggtitle("Basis2")) %>% print()
          })

#' @name plotting-themes
#' @aliases LinePlotTheme
#' @aliases EmptyTheme
#' @title Plotting themes
#' @description Formats a ggplot object for neat plotting.
#' @return Object of class \code{ggplot}
#' @export
#' @details \code{LinePlotTheme()} creates \code{ggplot} object with a white background, a relatively large font, and grid lines. \code{EmptyTheme()} on the other hand creates a \code{ggplot} object with no axes or legends.
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
    return(g)
}

#' @rdname plotting-themes
#' @export
EmptyTheme <- function() {
    g <- ggplot() +  theme(panel.background = element_rect(fill='white', colour='white'),
                           panel.grid=element_blank(),axis.ticks=element_blank(),
                           panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                           axis.text.x=element_blank(),axis.text.y=element_blank())
    return (g)
}

#####################################################
############# NOT EXPORTED ##########################
#####################################################

## Returns points ona circle with a given centre and diameter (for plotting it)
.circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2                            # radius
    tt <- seq(0,2*pi,length.out = npoints)      # default 100 points on circle
    xx <- center[1] + r * cos(tt)               # x values
    yy <- center[2] + r * sin(tt)               # y values
    return(data.frame(x = xx, y = yy))          # return in data frame
}

## This function ensures that there are no lines crossing the map due to countries traversing the
## -180, +180 boundary
.homogenise_maps <- function(worldmap) {
    group <- long <- prob <- NULL # suppress bindings note

    ## Take the world map, see which countries are "problematic" (prob == 1)
    ## and consider only those countries
    W <-worldmap %>%
        group_by(group) %>%
        summarise(prob = (max(long) > 180 | min(long) < -180)) %>%
        filter(prob==1)

    ## For each problematic country
    for(i in W$group) {
        this_country <- filter(worldmap,group == i)              # subset this country from worldmap
        CA <- filter(this_country, long >= 180 | long <= -180)   # find the problematic coordinates
        CB <- filter(this_country, long < 180 & long > -180)     # find the "OK" coordinates
        CA$group <- CA$group + 10000                             # put the problematic coordinates into a new group
        CB$group <- CB$group + 10001                             # put the new coordinates into a new group

        if(max(CA$long) >= 180) {                                # Shift all problematic longitudes that are too
            CA$long <- CA$long - 360                             # large to be within the [-180,180] range
        } else if(min(CA$long) <= -180) {                        # Same but for longitudes that are too small
            CA$long <- CA$long + 360
        }
        CA <- CA %>% filter(abs(long) <= 179.99)                 # If there are still problematic longitudes
        # just remove them
        worldmap <- rbind(worldmap,CA,CB)                        # add these to the world map
    }
    worldmap <- filter(worldmap,!(group %in% W$group))           # remove the problematic countries
    return(worldmap)                                             # return fixed world map
}



## Plotting functions: These may be removed, haven't decided yet. 

#' Plot non-Gaussian data and predictions.
#'
#' @inheritParams .plot_map
#' @inheritParams .plot_data
#' @param zdf A \code{dataframe} containing spatial coordinates (named "x" and "y") and the value of the observations (named "z").
#' @param newdata A \code{dataframe} containing the sspatial coordinates (named "x" and "y"), and the predictions and uncertainty quantification.
#' @return A list of ggplot() objects.
.plot_all <- function(newdata = NULL, zdf = NULL, response = NULL, x = "x", y = "y") {
    
    plots <- list()
    
    ## data plot
    if (!is.null(zdf)) plots$data <- .plot_data(zdf, response = response) + labs(colour = 'Z')
    
    ## Prediction and uncertainty plots
    if (!is.null(newdata)) {
        
        if ("p_Y" %in% names(newdata)) plots$p_Y <- .plot_map(newdata, col = "p_Y", x = x, y = y) + labs(fill = expression(widehat(p)[Y]["|"][bold(Z)]))
        if ("RMSPE_Y" %in% names(newdata)) plots$RMSPE_Y <- .plot_map(newdata, x = x, y = y, col = "RMSPE_Y", uncertaintyMap = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[Y]["|"][bold(Z)], Y))))
        if (all(c("Y_percentile_05","Y_percentile_95") %in% names(newdata))) {
            newdata$interval_90 <- newdata$Y_percentile_95 - newdata$Y_percentile_05 
            plots$interval_90_Y <-  .plot_map(newdata, x = x, y = y, col = "interval_90", uncertaintyMap = TRUE) + labs(fill = expression("90% Central \nInterval Width:" ~ Y))
        }
        
        if ("p_mu" %in% names(newdata)) plots$p_mu <- .plot_map(newdata, x = x, y = y, col = "p_mu") + labs(fill = expression(widehat(p)[mu]["|"][bold(Z)]))
        if ("RMSPE_mu" %in% names(newdata)) plots$RMSPE_mu <- .plot_map(newdata, x = x, y = y, col = "RMSPE_mu", uncertaintyMap = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[mu]["|"][bold(Z)], mu))))
        if (all(c("mu_percentile_05","mu_percentile_95") %in% names(newdata))) {
            newdata$interval_90 <- newdata$mu_percentile_95 - newdata$mu_percentile_05 
            plots$interval_90_mu <-  .plot_map(newdata, x = x, y = y, col = "interval_90", uncertaintyMap = TRUE) + labs(fill = expression("90% Central \nInterval Width:" ~ mu))
        }
        
        if ("p_prob" %in% names(newdata)) plots$p_prob <- .plot_map(newdata, x = x, y = y, col = "p_prob", diverging = TRUE, midpoint = 0.5) + labs(fill = expression(widehat(p)[pi]["|"][bold(Z)]))
        if ("RMSPE_prob" %in% names(newdata)) plots$RMSPE_prob <- .plot_map(newdata, x = x, y = y, col = "RMSPE_prob", uncertaintyMap = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[pi]["|"][bold(Z)], pi))))
        if (all(c("prob_percentile_05","prob_percentile_95") %in% names(newdata))) {
            newdata$interval_90 <- newdata$prob_percentile_95 - newdata$prob_percentile_05 
            plots$interval_90_prob <-  .plot_map(newdata, x = x, y = y, col = "interval_90", uncertaintyMap = TRUE) + labs(fill = expression("90% Central \nInterval Width:" ~ pi))
        }    
        
        if ("p_Z" %in% names(newdata)) plots$p_Z <- .plot_map(newdata, x = x, y = y, col = "p_Z") + labs(fill = expression(widehat(p)[Z]["|"][bold(Z)])) 
        # if ("p_Z_empirical" %in% names(newdata)) plots$p_Z_empirical <- .plot_map(newdata, col = "p_Z_empirical") + labs(fill = expression(widehat(p)[Z]["|"][bold(Z)]~" (empirical)"))
        if ("RMSPE_Z" %in% names(newdata)) plots$RMSPE_Z <- .plot_map(newdata, x = x, y = y, col = "RMSPE_Z", uncertaintyMap = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[Z]["|"][bold(Z)], Z))))
        if (all(c("Z_percentile_05","Z_percentile_95") %in% names(newdata))) {
            newdata$interval_90 <- newdata$Z_percentile_95 - newdata$Z_percentile_05 
            plots$interval_90_Z <-  .plot_map(newdata, x = x, y = y, col = "interval_90", uncertaintyMap = TRUE) + labs(fill = expression("90% Central \nInterval Width:" ~ Z))
        }    
    }
    
    return(plots)
}

#' Plot 2D spatial data
#'
#' Plot 2D spatial data using a diverging red and blue palette if 
#' \code{response == "bernoulli"}, and a spectral palette otherwise.
#'
#' @param zdf A \code{dataframe} containing spatial coordinates and the value of the response variable.
#' @param response A character indicating the assumed response distribution.
#' @param point_size Size of plotted points.
#' @param lim Controls the limits of the colour scale.
#' @param z A \code{string} indicating the data column.
#' @param x A \code{string} indicating the name of the column containing x-locations.
#' @param y A \code{string} indicating the name of the column containing the y-locations.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{.plot_map}}, \code{\link{.plot_all}}
.plot_data <- function(zdf,
                       response,
                       point_size = 1,
                       lim = range(zdf[["z"]]), 
                       z = "z", x = "x", y = "y"){
    
    ## Make ggplot() object
    p <- ggplot(zdf) +
        geom_point(aes_string(x = x, y = y, colour = z), size = point_size) +
        theme_bw() + coord_fixed()
    
    ## Colour scale
    if(response == "bernoulli"){
        p <- p + scale_colour_gradient2(low = "blue", high = "red", midpoint = 0.5, breaks = c(0, 1), guide = "legend")
    } else{
        p <- p + scale_colour_distiller(palette="Spectral", limits = lim)
    }
    
    return(p)
}



#' Plot of spatial process.
#'
#' @inheritParams .plot_data
#' @param df A \code{dataframe} containing spatial coordinates (named "x" and "y")
#' and the value of the process (whose name is specified by the \code{col} argument).
#' @param col A \code{string} indicating the name of the column containing the process values.
#' @param diverging \code{logical}, indicating whether a diverging palette should be used.
#' @param low Low colour (only applicable if \code{diverging == TRUE}).
#' @param mid Mid colour (only applicable if \code{diverging == TRUE}).
#' @param high High colour (only applicable if \code{diverging == TRUE}).
#' @param midpoint Point at which colour \code{mid} is assigned. This is useful for
#' uncertainty ratio maps (in which case a good midpoint is 1) and probability maps
#' (in which case a good midpoint is 0.5).
#' @param uncertaintyMap Logical indicating whether to use uncertainty colour scale.
#' @return A \code{ggplot} object of the spatial process.
#' @seealso \code{\link{.plot_data}}, \code{\link{.plot_all}}
.plot_map <- function(df, col, x = "x", y = "y", 
                      lim = range(df[[col]]), 
                      diverging = FALSE,
                      low = "blue", mid = "white", high = "red", midpoint = mean(df[[col]]),
                      uncertaintyMap = FALSE){
    
    ## Basic ggplot() object
    p <- ggplot(df) + geom_tile(aes_string(x = x, y = y, fill = col)) +
        theme_bw() + coord_fixed()
    
    ## Colour scale
    if (uncertaintyMap == TRUE & diverging == TRUE) {
        p <- p + scale_fill_gradient2(low = "#01665E", mid = "white", midpoint = midpoint, high = "#8C510A", limits = lim)
    } else if (uncertaintyMap == TRUE) {
        p <- p + scale_fill_distiller(palette = "BrBG", direction = -1, limits = lim)
    } else if (diverging == TRUE) {
        p <- p + scale_fill_gradient2(low = low, mid = mid, midpoint = midpoint, high = high, limits = lim)
    } else{
        p <- p + scale_fill_distiller(palette="Spectral", limits = lim)
    }
    
    return(p)
}
