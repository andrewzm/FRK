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

#' Plot predictions from FRK analysis. 
#' 
#' @rdname SRE
#' @param x \code{SRE} object 
#' @param y result of calling \code{predict} on an \code{SRE} object 
#' @param zdf a \code{data.frame}, \code{SpatialPointsDataFrame}, or \code{SpatialPolygonsDataFrame} containing the observations
#' @param map_layer a \code{ggplot} layer to add below the plotted layer, often a map
#' @return A list of \code{ggplot} objects consisting of the observed data, predictions, and standard errors. This list can then be supplied to, for example, \code{ggpubr::ggarrange()}.
#' @export
setMethod("plot", signature(x = "SRE"), function(x, y, zdf = NULL, map_layer = NULL, ...) {
    
    SRE_model <- x
    pred_object <- y
    
    ## Check that pred_object is a result of a call to predict()
    if (SRE_model@method == "TMB") {
        if (!is(pred_object, "list")) 
            stop("Since method = 'TMB', y (the prediction object) should be a list")
        if (is.null(pred_object$newdata)) 
            stop("Since method = 'TMB', y (the prediction object) should be a list with an element called newdata")
        pred_object <- pred_object$newdata 
    } else if (SRE_model@method == "EM") {
        if(!is(pred_object, "Spatial") && !is(pred_object, "ST")) 
            stop("Since method = 'EM', y (the prediction object) should be a Spatial* or ST* object")
    }
        
    
    
    ## Plot the predictions and UQ 
    plots <- .plot_predictions_and_UQ(SRE_model, pred_object, map_layer, ...)
    
    
    
    ## Plot the data if it was provided
    if (!is.null(zdf)) {
        ## Extract name of response variable we wish to plot
        response_name <- all.vars(SRE_model@f)[1]
        plots$data <- .plot_data(zdf, response_name, map_layer, ...)
    }
        
    

    
    return(plots)
})

## Restructure without having to pass zdf to .plot_data. This is so I can use 
## .plot_data on any Spatial* object. 
.plot_predictions_and_UQ <- function(SRE_model, pred_object, map_layer = NULL, ...) {
    
    plots <- list()
    
    ## FIXME: the only thing we need here from SRE_model is method (if we can get time_name from pred_object)
    
    ## Get the coordinate names from the pred_object
    ## TODO: See if it makes sense to base the following segment on pred_object rather 
    ## than BAUs. 
    ## FIXME: there is a lot of repetition between this function and .plot_data.
    coord_names <- coordnames(SRE_model@BAUs)
    ## Remove the coordinate corresponding to time. coord_names is just spatial.
    if(is(SRE_model@BAUs, "ST")) {
        time_name <- coord_names[3]
        coord_names <- coord_names[1:2]
    } else {
        time_name <- NULL
    }
    
    ## Plot the predictions and uncertainty.
    ## Convert to data.frame for plotting. Do this here so we only have to do it 
    ## once, as this can be expensive for large objects. 
    if (is(pred_object, "SpatialPixelsDataFrame")) {
        df <- data.frame(pred_object)
        sp_type <- "pixels"
    } else if (is(pred_object, "SpatialPolygonsDataFrame")) {
        df <- .SpatialPolygonsDataFrame_to_df(pred_object)
        sp_type <- "polygons"
    } else if (is(pred_object, "STFDF")) {
        if (is(pred_object@sp, "SpatialPolygons")) {
            sp_type <- "polygons"
        } else if (is(pred_object@sp, "SpatialPixels")) {
            sp_type <- "pixels"
        }
        df <- STFDF_to_df(pred_object) 
    } else {
        stop("Class of y (object from a call to predict) is not recognised by plot().")
    }
    
    # df, col, coord_names, time_name, sp_type, map_layer
    
    if (SRE_model@method == "EM") {
        plots$mu <- .plot(df, "mu", coord_names, time_name, sp_type, map_layer, ...) 
        plots$se <- .plot(df, "sd", coord_names, time_name, sp_type, map_layer, uncertainty = TRUE, ...) + labs(fill = "se")
        
    } else if (SRE_model@method == "TMB") {
        
        ## Unfortunately, we don't know the value of type specified in predict().
        ## We will just check the existence of quantities using if statements. 
        ## Usually, people want prediction interval widths, but the output of predict
        ## doesn't give interval widths directly. Instead, it gives user-specified
        ## percentiles, so here we construct some interval widths assuming the user
        ## has specified the 5th and 95th percentiles (default).
        
        if ("p_Y" %in% names(df)) plots$p_Y <- .plot(df, "p_Y", coord_names, time_name, sp_type, map_layer, ...) + labs(fill = expression(widehat(p)[Y]["|"][bold(Z)]))
        if ("RMSPE_Y" %in% names(df)) plots$RMSPE_Y <- .plot(df, "RMSPE_Y", coord_names, time_name, sp_type, map_layer, uncertainty = TRUE, map_layer = map_layer, ...) + labs(fill = expression(sqrt(MSPE(widehat(p)[Y]["|"][bold(Z)], Y))))
        if (all(c("Y_percentile_5","Y_percentile_95") %in% names(df))) {
            df$interval_90 <- df$Y_percentile_95 - df$Y_percentile_5 
            plots$interval_90_Y <-  .plot(df, "interval_90", coord_names, time_name, sp_type, map_layer, uncertainty = TRUE, ...) + 
                labs(fill = eval(bquote(expression(
                    "Width of 90%\npredictive interval\nfor latent process Y(" *"\U00B7)"
                )))) 
        }
        
        if ("p_mu" %in% names(df)) plots$p_mu <- .plot(df, "p_mu", coord_names, time_name, sp_type, map_layer, ...) + labs(fill = expression(widehat(p)[mu]["|"][bold(Z)]))
        if ("RMSPE_mu" %in% names(df)) plots$RMSPE_mu <- .plot(df, "RMSPE_mu", coord_names, time_name, sp_type, map_layer, uncertainty = TRUE, ...) + labs(fill = expression(sqrt(MSPE(widehat(p)[mu]["|"][bold(Z)], mu))))
        if (all(c("mu_percentile_5","mu_percentile_95") %in% names(df))) {
            df$interval_90 <- df$mu_percentile_95 - df$mu_percentile_5 
            plots$interval_90_mu <-  .plot(df, "interval_90", coord_names, time_name, sp_type, map_layer, uncertainty = TRUE, ...) + 
                labs(fill = eval(bquote(expression(
                    "Width of 90%\npredictive interval\nfor mean process " *mu *"(\U00B7)"
                )))) 
        }
        
        if ("p_prob" %in% names(df)) plots$p_prob <- .plot(df, "p_prob", coord_names, time_name, sp_type, map_layer, ...) + labs(fill = expression(widehat(p)[pi]["|"][bold(Z)]))
        if ("RMSPE_prob" %in% names(df)) plots$RMSPE_prob <- .plot(df, "RMSPE_prob", coord_names, time_name, sp_type, map_layer, uncertainty = TRUE, ...) + labs(fill = expression(sqrt(MSPE(widehat(p)[pi]["|"][bold(Z)], pi))))
        if (all(c("prob_percentile_5","prob_percentile_95") %in% names(df))) {
            df$interval_90 <- df$prob_percentile_95 - df$prob_percentile_5 
            plots$interval_90_prob <-  .plot(df, "interval_90", coord_names, time_name, sp_type, map_layer, uncertainty = TRUE, ...) + 
                labs(fill = eval(bquote(expression(
                    "Width of 90%\npredictive interval\nfor probability process " *pi *"(\U00B7)"
                )))) 
        }   
        
        if ("p_Z" %in% names(df)) plots$p_Z <- .plot(df, "p_Z", coord_names, time_name, sp_type, map_layer, ...) + labs(fill = expression(widehat(p)[Z]["|"][bold(Z)]))
        if ("RMSPE_Z" %in% names(df)) plots$RMSPE_Z <- .plot(df, "RMSPE_Z", coord_names, time_name, sp_type, map_layer, uncertainty = TRUE, ...) + labs(fill = expression(sqrt(MSPE(widehat(p)[Z]["|"][bold(Z)], Z))))
        if (all(c("Z_percentile_5","Z_percentile_95") %in% names(df))) {
            df$interval_90 <- df$Z_percentile_95 - df$Z_percentile_5 
            plots$interval_90_Z <-  .plot(df, "interval_90", coord_names, time_name, sp_type, map_layer, uncertainty = TRUE, ...) + 
                labs(fill = eval(bquote(expression(
                    "Width of 90%\npredictive interval\nfor data process " *Z *"(\U00B7)"
                )))) 
        }
    }
    
    return(plots)
    
}


## Plot of predictions or uncertainty quantification.
.plot <- function(df, col, coord_names, time_name, sp_type, map_layer, uncertainty = FALSE, ...){

    ## Remove duplicate columns (can sometimes happen when we convert the Spatial* 
    ## object, if the coordinates are already present)
    df <- df[, unique(colnames(df))]
    
    if (length(coord_names) != 2) 
        stop("plot() only implemented for the 2D space")
    
    ## Edit the time column so that the facets display t = ...
    if (!is.null(time_name)) 
        df[, time_name] <- factor(
            df[, time_name], ordered = TRUE,
            labels = paste(time_name, sort(unique(df[, time_name])), sep = " = ")
        )
        
    ## Basic plot
    gg <- ggplot(data = df, aes_string(x = coord_names[1], y = coord_names[2])) + 
        theme_bw() + coord_fixed()
    
    ## Add map_layer (if present)
    if (!is.null(map_layer)) 
        gg <- gg + map_layer
    
    ## Plot based on data type
    if (sp_type == "points") { # this is only for the observations
        gg <- gg + geom_point(aes_string(colour = col), ...) + scale_colour_distiller(palette = "Spectral")
    } else if (sp_type == "pixels") {
        gg <- gg + geom_raster(aes_string(fill = col), ...)
    } else if (sp_type == "polygons") {
        gg <- gg + geom_polygon(aes_string(group = "id", fill = col), ...) 
    } 
    
    if (!is.null(time_name))
        gg <- gg + facet_wrap(as.formula(paste("~", time_name)))
    
    ## fill scale
    if (uncertainty) {
        gg <- gg + scale_fill_distiller(palette = "BrBG", direction = -1)
    } else {
        gg <- gg + scale_fill_distiller(palette = "Spectral")
    }
    
    return(gg)
}





## TODO: could make this an exported function for plotting data. 
## Would be useful in the Sydney example for making the SA2 plots. 
## We would just need an additional argument, col, which allows one to specify 
## which column to plot. 
.plot_data <- function(zdf, response_name, map_layer = NULL, ...){
    
    if (!is(zdf, "Spatial") && !is(zdf, "STFDF")) 
        stop("zdf should be a Spatial*DataFrame or STFDF")
    
    ## Get the coordinate names 
    coord_names <- coordnames(zdf)
    
    ## Remove the coordinate corresponding to time. coord_names is just spatial.
    if(is(zdf, "ST")) {
        time_name <- coord_names[3]
        coord_names <- coord_names[1:2]
    } else {
        time_name <- NULL
    }
    
    if (is(zdf, "SpatialPointsDataFrame")) {
        zdf <- data.frame(zdf)
        sp_type <- "points"
    } else if (is(zdf, "SpatialPolygonsDataFrame")) {
        zdf <- .SpatialPolygonsDataFrame_to_df(zdf)
        sp_type <- "polygons"
    } else if (is(zdf, "STIDF")) {
        stop("Plotting of STIDF not yet implemented.")
        zdf <- data.frame(zdf)
        sp_type <- "points"
    } else if (is(zdf, "STFDF")) {
        if (is(zdf@sp, "SpatialPolygons")) {
            sp_type <- "polygons"
        } else if (is(zdf@sp, "SpatialPixels")) {
            sp_type <- "pixels"
        }
        zdf <- STFDF_to_df(zdf) 
    } else {
        stop("Class of zdf is not recognised by plot().")
    }
    
    return(.plot(zdf, response_name, coord_names, time_name, sp_type, map_layer, ...))
}