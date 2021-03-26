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

# #' Plot predictions from FRK analysis.
# #' 
# #' @rdname SRE
# #' @param object \code{SRE} object
# #' @param Pred result of a calling \code{predict} on an \code{SRE} object
# #' @param zdf a \code{data.frame}, \code{SpatialPointsDataFrame}, or \code{SpatialPolygonsDataFrame} containing the observations
# #' @return A list of \code{ggplot} objects consisting of the observed data, predictions, and standard errors. This list can then be supplied to, for example, \code{ggpubr::ggarrange()}.
# #' @export
# # SRE.plot <- function(object, Pred, zdf = NULL) {
# setMethod("plot", signature(object = "SRE"), function(object, Pred, zdf = NULL) {
#     
#     plots <- list() # initialise plot list
#     
#     if (object@method == "TMB")
#         Pred <- Pred$newdata # If method = TMB, Pred is a list
#     
#     if(is(Pred, "ST"))
#         stop("SRE.plot only implemented for spatial applications")
#     
#     ## Extract names of coordinates
#     coord_names <- coordnames(Pred) 
#     
#     if (length(coord_names) != 2) 
#         stop("SRE.plot only implemented for two-dimensional Euclidean space")
#     
#     if (!is.null(zdf)) 
#         plots$data <- .plot_data(zdf, object, coord_names)
#     
#     if (is(Pred, "SpatialPixelsDataFrame")) {
#         xy <- data.frame(Pred)
#         sp_type <- "pixels"
#     } else if (is(Pred, "SpatialPolygonsDataFrame")) {
#         xy <- .SpatialPolygonsDataFrame_to_df(Pred)
#         sp_type <- "polygons"
#     } else {
#         stop("Class of Pred is not recognised by SRE.plot().")
#     }
# 
#     if (object@method == "EM") {
#         plots$mu <- .plot_map(xy, coord_names, col = "mu", sp_type)
#         plots$se <- .plot_map(xy, coord_names, col = "sd", sp_type, uncertainty = TRUE)
#         
#     } else if (object@method == "TMB") {
#         
#         ## Unfortunately, we don't know the value of type specified in predict().
#         ## We will just check the existence of quantities using if statements. 
#         ## Usually, people want prediction interval widths, but the output of predict
#         ## doesn't give interval widths directly. Instead, it gives user-specified
#         ## percentiles, so here we construct some interval widths assuming the user
#         ## has specified the 5th and 95th percentiles (default).
#         
#         if ("p_Y" %in% names(xy)) plots$p_Y <- .plot_map(xy, coord_names, col = "p_Y", sp_type) + labs(fill = expression(widehat(p)[Y]["|"][bold(Z)]))
#         if ("RMSPE_Y" %in% names(xy)) plots$RMSPE_Y <- .plot_map(xy, coord_names, col = "RMSPE_Y", sp_type, uncertainty = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[Y]["|"][bold(Z)], Y))))
#         if (all(c("Y_percentile_5","Y_percentile_95") %in% names(xy))) {
#             xy$interval_90 <- xy$Y_percentile_95 - xy$Y_percentile_5 
#             plots$interval_90_Y <-  .plot_map(xy, coord_names, col = "interval_90", sp_type, uncertainty = TRUE) + 
#                 labs(fill = eval(bquote(expression(
#                     "Width of 90%\npredictive interval\nfor latent process Y(" *"\U00B7)"
#                 )))) 
#         }
#         
#         if ("p_mu" %in% names(xy)) plots$p_mu <- .plot_map(xy, coord_names, col = "p_mu", sp_type) + labs(fill = expression(widehat(p)[mu]["|"][bold(Z)]))
#         if ("RMSPE_mu" %in% names(xy)) plots$RMSPE_mu <- .plot_map(xy, coord_names, col = "RMSPE_mu", sp_type, uncertainty = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[mu]["|"][bold(Z)], mu))))
#         if (all(c("mu_percentile_5","mu_percentile_95") %in% names(xy))) {
#             xy$interval_90 <- xy$mu_percentile_95 - xy$mu_percentile_5 
#             plots$interval_90_mu <-  .plot_map(xy, coord_names, col = "interval_90", sp_type, uncertainty = TRUE) + 
#                 labs(fill = eval(bquote(expression(
#                     "Width of 90%\npredictive interval\nfor mean process " *mu *"(\U00B7)"
#                 )))) 
#         }
#         
#         if ("p_prob" %in% names(xy)) plots$p_prob <- .plot_map(xy, coord_names, col = "p_prob", sp_type) + labs(fill = expression(widehat(p)[pi]["|"][bold(Z)]))
#         if ("RMSPE_prob" %in% names(xy)) plots$RMSPE_prob <- .plot_map(xy, coord_names, col = "RMSPE_prob", sp_type, uncertainty = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[pi]["|"][bold(Z)], pi))))
#         if (all(c("prob_percentile_5","prob_percentile_95") %in% names(xy))) {
#             xy$interval_90 <- xy$prob_percentile_95 - xy$prob_percentile_5 
#             plots$interval_90_prob <-  .plot_map(xy, coord_names, col = "interval_90", sp_type, uncertainty = TRUE) + 
#                 labs(fill = eval(bquote(expression(
#                     "Width of 90%\npredictive interval\nfor probability process " *pi *"(\U00B7)"
#                 )))) 
#         }   
#         
#         if ("p_Z" %in% names(xy)) plots$p_Z <- .plot_map(xy, coord_names, col = "p_Z", sp_type) + labs(fill = expression(widehat(p)[Z]["|"][bold(Z)]))
#         if ("RMSPE_Z" %in% names(xy)) plots$RMSPE_Z <- .plot_map(xy, coord_names, col = "RMSPE_Z", sp_type, uncertainty = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[Z]["|"][bold(Z)], Z))))
#         if (all(c("Z_percentile_5","Z_percentile_95") %in% names(xy))) {
#             xy$interval_90 <- xy$Z_percentile_95 - xy$Z_percentile_5 
#             plots$interval_90_Z <-  .plot_map(xy, coord_names, col = "interval_90", sp_type, uncertainty = TRUE) + 
#                 labs(fill = eval(bquote(expression(
#                     "Width of 90%\npredictive interval\nfor data process " *Z *"(\U00B7)"
#                 )))) 
#         }
#     }
#     
#     return(plots)
# })


#' Plot predictions from FRK analysis. 
#' 
#' @rdname SRE
#' @param object \code{SRE} object 
#' @param y result of calling \code{predict} on an \code{SRE} object 
#' @param zdf a \code{data.frame}, \code{SpatialPointsDataFrame}, or \code{SpatialPolygonsDataFrame} containing the observations
#' @return A list of \code{ggplot} objects consisting of the observed data, predictions, and standard errors. This list can then be supplied to, for example, \code{ggpubr::ggarrange()}.
#' @export
setMethod("plot", signature(x = "SRE"), function(x, y, zdf = NULL) {
    
    plots <- list() # initialise plot list
    
    if (x@method == "TMB")
        y <- y$newdata # If method = TMB, y is a list
    
    if(is(y, "ST"))
        stop("plot only implemented for spatial applications")
    
    ## Extract names of coordinates
    coord_names <- coordnames(y) 
    
    if (length(coord_names) != 2) 
        stop("plot only implemented for two-dimensional Euclidean space")
    
    if (!is.null(zdf)) 
        plots$data <- .plot_data(zdf, x, coord_names)
    
    if (is(y, "SpatialPixelsDataFrame")) {
        df <- data.frame(y)
        sp_type <- "pixels"
    } else if (is(y, "SpatialPolygonsDataFrame")) {
        df <- .SpatialPolygonsDataFrame_to_df(y)
        sp_type <- "polygons"
    } else {
        stop("Class of y is not recognised by SRE.plot().")
    }
    
    if (x@method == "EM") {
        plots$mu <- .plot_map(df, coord_names, col = "mu", sp_type)
        plots$se <- .plot_map(df, coord_names, col = "sd", sp_type, uncertainty = TRUE)
        
    } else if (x@method == "TMB") {
        
        ## Unfortunately, we don't know the value of type specified in predict().
        ## We will just check the existence of quantities using if statements. 
        ## Usually, people want prediction interval widths, but the output of predict
        ## doesn't give interval widths directly. Instead, it gives user-specified
        ## percentiles, so here we construct some interval widths assuming the user
        ## has specified the 5th and 95th percentiles (default).
        
        if ("p_Y" %in% names(df)) plots$p_Y <- .plot_map(df, coord_names, col = "p_Y", sp_type) + labs(fill = expression(widehat(p)[Y]["|"][bold(Z)]))
        if ("RMSPE_Y" %in% names(df)) plots$RMSPE_Y <- .plot_map(df, coord_names, col = "RMSPE_Y", sp_type, uncertainty = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[Y]["|"][bold(Z)], Y))))
        if (all(c("Y_percentile_5","Y_percentile_95") %in% names(df))) {
            df$interval_90 <- df$Y_percentile_95 - df$Y_percentile_5 
            plots$interval_90_Y <-  .plot_map(df, coord_names, col = "interval_90", sp_type, uncertainty = TRUE) + 
                labs(fill = eval(bquote(expression(
                    "Width of 90%\npredictive interval\nfor latent process Y(" *"\U00B7)"
                )))) 
        }
        
        if ("p_mu" %in% names(df)) plots$p_mu <- .plot_map(df, coord_names, col = "p_mu", sp_type) + labs(fill = expression(widehat(p)[mu]["|"][bold(Z)]))
        if ("RMSPE_mu" %in% names(df)) plots$RMSPE_mu <- .plot_map(df, coord_names, col = "RMSPE_mu", sp_type, uncertainty = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[mu]["|"][bold(Z)], mu))))
        if (all(c("mu_percentile_5","mu_percentile_95") %in% names(df))) {
            df$interval_90 <- df$mu_percentile_95 - df$mu_percentile_5 
            plots$interval_90_mu <-  .plot_map(df, coord_names, col = "interval_90", sp_type, uncertainty = TRUE) + 
                labs(fill = eval(bquote(expression(
                    "Width of 90%\npredictive interval\nfor mean process " *mu *"(\U00B7)"
                )))) 
        }
        
        if ("p_prob" %in% names(df)) plots$p_prob <- .plot_map(df, coord_names, col = "p_prob", sp_type) + labs(fill = expression(widehat(p)[pi]["|"][bold(Z)]))
        if ("RMSPE_prob" %in% names(df)) plots$RMSPE_prob <- .plot_map(df, coord_names, col = "RMSPE_prob", sp_type, uncertainty = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[pi]["|"][bold(Z)], pi))))
        if (all(c("prob_percentile_5","prob_percentile_95") %in% names(df))) {
            df$interval_90 <- df$prob_percentile_95 - df$prob_percentile_5 
            plots$interval_90_prob <-  .plot_map(df, coord_names, col = "interval_90", sp_type, uncertainty = TRUE) + 
                labs(fill = eval(bquote(expression(
                    "Width of 90%\npredictive interval\nfor probability process " *pi *"(\U00B7)"
                )))) 
        }   
        
        if ("p_Z" %in% names(df)) plots$p_Z <- .plot_map(df, coord_names, col = "p_Z", sp_type) + labs(fill = expression(widehat(p)[Z]["|"][bold(Z)]))
        if ("RMSPE_Z" %in% names(df)) plots$RMSPE_Z <- .plot_map(df, coord_names, col = "RMSPE_Z", sp_type, uncertainty = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[Z]["|"][bold(Z)], Z))))
        if (all(c("Z_percentile_5","Z_percentile_95") %in% names(df))) {
            df$interval_90 <- df$Z_percentile_95 - df$Z_percentile_5 
            plots$interval_90_Z <-  .plot_map(df, coord_names, col = "interval_90", sp_type, uncertainty = TRUE) + 
                labs(fill = eval(bquote(expression(
                    "Width of 90%\npredictive interval\nfor data process " *Z *"(\U00B7)"
                )))) 
        }
    }
    
    return(plots)
})



## Plot of predictions or uncertainty quantification.
.plot_map <- function(df, coord_names, col, sp_type, uncertainty = FALSE){
    
    ## Basic plot
    gg <- ggplot(data = df, aes_string(x = coord_names[1], y = coord_names[2])) + 
        theme_bw() + coord_fixed()
    
    ## Plot based on data type
    if (sp_type == "pixels") {
        gg <- gg + geom_raster(aes_string(fill = col))
    } else if (sp_type == "polygons") {
        gg <- gg + geom_polygon(aes_string(group = "id", fill = col), 
                                colour = "black"
                                ) 
    }
    
    ## Colour scale
    if (uncertainty) {
        gg <- gg + scale_fill_distiller(palette = "BrBG", direction = -1)
    } else {
        gg <- gg + scale_fill_distiller(palette = "Spectral")
    }
    
    
    return(gg)
}


## Plot 2D spatial data
.plot_data <- function(zdf, object, coord_names){
    
    ## Extract name of response variable so we can plot 
    response_name <- all.vars(object@f)[1]
    
    if (is(zdf, "data.frame") || is(zdf, "SpatialPointsDataFrame")) {
        
        zdf <- data.frame(zdf)
        data_plot <- ggplot(zdf) + geom_point(
            aes_string(x = coord_names[1], y = coord_names[2], colour = response_name)) + 
            scale_colour_distiller(palette = "Spectral") + theme_bw() + coord_fixed()
        
    } else if (is(zdf, "SpatialPolygonsDataFrame")) {
        
        zdf <- .SpatialPolygonsDataFrame_to_df(zdf)
        data_plot <- ggplot(zdf) +
            geom_polygon(aes_string(
                coord_names[1], y = coord_names[2], group = "id", fill = response_name), 
                colour = "black") +
            scale_fill_distiller(palette = "Spectral") +
            coord_fixed() + theme_bw()
        
    } else {
        stop("Class of zdf is not recognised by SRE.plot().")
    }
    
    return(data_plot)
}