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
#' @details This function uses \code{ggplot2::map_data()} in order to create a world map. Since, by default, this creates lines crossing the world at the (-180,180) longitude boundary, the function \code{.homogenise_maps()} is used to split the polygons at this boundary into two. If \code{inc_border} is TRUE, then a border is drawn around the lon-lat space; this option is most useful for projections that do not yield rectangular plots (e.g., the sinusoidal global projection).
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
    g + geom_path(data = worldmap, aes(x=long, y=lat, group=group), colour="black",linewidth=0.1)
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
    g <- ggplot() + .EmptyTheme_theme_only()
    return (g)
}

## this is so I can add it to already-constructed plots
.EmptyTheme_theme_only <- function() {
    theme(panel.background = element_rect(fill='white', colour='white'),
          panel.grid=element_blank(),axis.ticks=element_blank(),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          axis.text.x=element_blank(),axis.text.y=element_blank())
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

# ---- FUNCTIONS SPECIFIC TO PLOT() -----

#' @rdname plot
#' @export 
setMethod("plot", signature(x = "SRE", y = "list"), function(x, y, ...) {
  
  if(!("newdata" %in% names(y))) stop("y should contain an element named 'newdata'")
  
  plot(x, y, ...)
})

#' @rdname plot
#' @export 
setMethod("plot", signature(x = "SRE", y = "STFDF"), function(x, y, ...) .plot_common(x, y, ...))

#' @rdname plot
#' @export 
setMethod("plot", signature(x = "SRE", y = "SpatialPointsDataFrame"), function(x, y, ...) .plot_common(x, y, ...))

#' @rdname plot
#' @export 
setMethod("plot", signature(x = "SRE", y = "SpatialPixelsDataFrame"), function(x, y, ...) .plot_common(x, y, ...))

#' @rdname plot
#' @export 
setMethod("plot", signature(x = "SRE", y = "SpatialPolygonsDataFrame"), function(x, y, ...) .plot_common(x, y, ...))


.plot_common <- function(object, newdata, ...) {
  
  ## Determine columns we need to plot, and which palette to use
  if (object@method == "EM") {
    
    ## Here we will just rename the columns produced by the EM side to 
    ## align with those of the TMB side.
    newdata$p_mu <- newdata$mu 
    newdata$RMSPE_mu <- newdata$sd
    column_names <- c("p_mu", "RMSPE_mu")
    palette <- c("Spectral", "BrBG")
    
  } else if (object@method == "TMB") {
    
    ## Unfortunately, we don't know the value of type specified in predict().
    ## We will just check the existence of quantities using if statements. 
    ## Usually, people want prediction interval widths, but the output of predict
    ## doesn't give interval widths directly. Instead, it gives user-specified
    ## percentiles, so here we construct some interval widths assuming the user
    ## has specified the 5th and 95th percentiles (default).
    
    ## If the 5th and 95th percentiles are present, construct the 90% predictive interval width. 
    if (all(c("Y_percentile_5","Y_percentile_95") %in% names(newdata@data))) 
      newdata@data$interval90_Y <- newdata$Y_percentile_95 - newdata$Y_percentile_5 
    if (all(c("mu_percentile_5","mu_percentile_95") %in% names(newdata@data))) 
      newdata@data$interval90_mu <- newdata@data$mu_percentile_95 - newdata@data$mu_percentile_5 
    if (all(c("prob_percentile_5","prob_percentile_95") %in% names(newdata@data))) 
      newdata@data$interval90_prob <- newdata@data$prob_percentile_95 - newdata@data$prob_percentile_5 
    if (all(c("Z_percentile_5","Z_percentile_95") %in% names(newdata@data))) 
      newdata@data$interval90_Z <- newdata@data$Z_percentile_95 - newdata@data$Z_percentile_5 
    
    ## Determine which quantities we are actually plotting
    QOI <- c("Y", "mu", "prob", "Z") # Possible Quantities of interest
    column_names <- paste0(rep(c("p_", "RMSPE_", "interval90_"), each = length(QOI)), QOI)
    palette <- rep(c("Spectral", "BrBG", "BrBG"), each = length(QOI)) 
    keep <- column_names %in% names(newdata@data)
    column_names <- column_names[keep]
    palette  <-  palette[keep] 
  }
  
  ## plot() decides the palette for each graphic based on whether it is 
  ## the term is related to prediction or prediction uncertainty quantification.
  ## This means that the palette argument cannot be provided by the user via 
  ## ... argument. An exception is if the user has specified palette = "nasa", 
  ## in which case we will use it for the prediction plots.
  user_args <- list(...)
  if (!is.null(user_args$palette) && user_args$palette == "nasa") {
    palette <- gsub(pattern = "Spectral", replacement = "nasa", palette)
  }
  
  ## Need to use do.call(), as it allows arguments to be provided as a list. 
  ## Otherwise, there is no way to remove the original value of palette supplied
  ## by the user. Remove potential clashes:
  user_args$palette      <- NULL
  user_args$column_names <- NULL
  args_list <- c(list(newdata = newdata, 
                      column_names = column_names, 
                      palette = palette), 
                 user_args)
  plots <- do.call(plot_spatial_or_ST, args_list)
  

  ## First we determine if the user has predicted over the BAUs or over arbitrary
  ## prediction regions. To assess this, we will test if newdata has the same 
  ## number of elements as the BAUs (this is not foolproof, but good enough).
  pred_over_BAUs <- length(object@BAUs) == length(newdata)
  
  ## Edit labels
  split_column_names <- strsplit(column_names, "_")
  names(split_column_names) <- column_names
  for (i in column_names) {
    plots[[i]] <- plots[[i]] + .custom_lab(split_column_names[[i]], pred_over_BAUs)
  }
  
  ## Now plot the data. Note that here we use a simple call to plot_spatial_or_ST(), 
  ## rather than via do.call(), because we have not constructed a custom palette 
  ## within this function.
  response_name <- all.vars(object@f)[1]
  if (is.data.frame(object@data)) {
    ## This will never happen, but add a check just in case
    warning("Cannot plot the data stored in object@data, as it is a data.frame")
  } else if (is.list(object@data)) {
    if (any(sapply(object@data, function(d) is(d, "STIDF")))) {
      if (length(object@data) == 1) {
        ## Plot the binned data over the BAUs. This is a bit difficult if any 
        ## observations are associated with multiple BAUs; don't do it in this case.
        Cmat_dgT <- .as(object@Cmat, "dgTMatrix")
        if (!all(tabulate(Cmat_dgT@i + 1) == 1) || any(tabulate(Cmat_dgT@j + 1) > 1)) {
         cat("The data set stored in object@data is of class STIDF. In this case, 
              we normally use the binned data in object@Z to plot over the BAUs. However,
              some observations are associated with multiple BAUs, or there 
              are some BAUs associated with multiple observations (probably because 
              average_in_BAU = FALSE). This complicates matters, so we will not plot the data.")
        } else {
          Z <- binned_data(object) 
          ## If NAs are present (i.e., some of the BAUs are unobserved), 
          ## tell the user. This will illuminate a warning thrown later as well.
          if (any(is.na(Z))) {
            cat("To plot the STIDF data provided in the SRE object, we use the binned data in object@Z to plot over the BAUs. The unobserved BAUs (i.e., those that are not associated with any elements of object@Z) are given a value of NA.\n")
          }
          object@BAUs@data[, response_name] <- Z
          data_plots <- plot_spatial_or_ST(object@BAUs, response_name, ...)
        }
      } else {
        ## We may have a combination of STIDF and STFDF, and we don't know which 
        ## dataset each element of object@Z is associated with. This is a bit 
        ## complicated, so we don't do it. 
        cat("Multiple data sets were used in the analysis, and because at 
                least one is of class STIDF, we will not plot the data (we don't have a method for this).\n")
      }
    } else {
      data_plots <- sapply(object@data, plot_spatial_or_ST, response_name, ...)
    }
  } else if (is(object@data, "SpatialPointsDataFrame") ||
             is(object@data, "SpatialPixelsDataFrame") ||
             is(object@data, "SpatialPolygonsDataFrame") ||
             is(object@data, "STFDF")) { 
    ## this shouldn't happen, but add just in case
    data_plots <- plot_spatial_or_ST(object@data, response_name, ...)
  } else {
    warning("Couldn't plot the data stored in object@data as the type was unrecognised.")
  }
  
  ## Give the data plots meaningful labels, and combine them with the prediction 
  ## and uncertainty quantification plots.
  if (exists("data_plots")) {
    if (length(data_plots) == 1) {
      names(data_plots) <- response_name
    } else {
      names(data_plots) <- paste(response_name, "dataset", 1:length(data_plots), sep = "_")
      data_legend_name <- paste(response_name, "\ndataset", 1:length(data_plots), sep = " ")
      for (i in 1:length(data_plots)) {
        data_plots[[i]] <- data_plots[[i]] + labs(colour = data_legend_name[i], fill = data_legend_name[i])
      }
    }
    plots <- c(data_plots, plots) 
  }
  
  return(plots)
}



.custom_lab <- function(v, pred_over_BAUs) {
    # v is a vector, with first entry indicating the type of plot we are making a 
    ## label for, and the second entry indicating the quantity of interest    
    type <- v[1]; x <- v[2]
    
    ## Set the unicode character for the quantity of interest
    ## See https://en.wikipedia.org/wiki/List_of_Unicode_characters#Greek_and_Coptic
    # for the a list of unicode characters.
    unicode <- if (x == "Y") "Y" else if (x == "mu") "\U03BC" else if (x == "prob") "\U03C0" else if (x == "Z") "Z*"

    if (pred_over_BAUs) {
      process <-  bquote(bold(.(unicode)))
    } else {
      process <-  bquote(bold(.(unicode))[P])
    }
    
    expectation <- bquote(paste("E(", .(process), " | ", bold(Z), ", ", bold("\U03B8"), ")"))
    
    ## Construct the labels
    ## NB: Added white space to ensure no overlap between label and the fill box 
    label <- if (type == "p") {
      # bquote(atop("Prediction    ", .(expectation))) # TODO add the following whitespace to expectation "    "
      top    <- "Prediction    "
      bottom <- bquote(paste(.(expectation) * "    "))
      bquote(atop(.(top), .(bottom))) 
    } else if (type == "RMSPE") {
        bquote(paste("RMSPE(" * .(expectation) * ", " * .(process),")        "))
    } else if (type == "interval90") {
        top    <- "90% prediction-    "
        bottom <- bquote(paste("interval width for " * .(process), "  "))
        bquote(atop(.(top), .(bottom)))
    }
    
    return(labs(fill = label))
}


# ---- MAIN PLOTTING FUNCTIONS ----


#' @rdname plot_spatial_or_ST
#' @export
setMethod("plot_spatial_or_ST", signature(newdata = "STFDF"), 
          function(newdata, column_names,  map_layer = NULL, 
                   subset_time = NULL, palette = "Spectral", 
                   plot_over_world = FALSE, labels_from_coordnames = TRUE, ...) {
      
      ## Get the coordinate names 
      coord_names <- coordnames(newdata)
      
      ## Remove the coordinate corresponding to time. coord_names is just spatial.
      time_name <- coord_names[3]
      coord_names <- coord_names[1:2]
      
      if (!is.null(subset_time)) 
          newdata <- newdata[, subset_time]
                  
      plots <- .plot_spatial_or_ST_common(
          newdata = newdata, column_names = column_names, coord_names = coord_names, 
          time_name = time_name, map_layer = map_layer, 
          palette = palette, plot_over_world = plot_over_world, 
          labels_from_coordnames = labels_from_coordnames,
          ...
          )
      return(plots)
})        


#' @rdname plot_spatial_or_ST
#' @export
setMethod("plot_spatial_or_ST", signature(newdata = "SpatialPointsDataFrame"), 
          function(newdata, column_names,  map_layer = NULL, 
                   subset_time = NULL, palette = "Spectral", 
                   plot_over_world = FALSE, labels_from_coordnames = TRUE, ...) {
    
    ## Get the coordinate names, and time_name is NULL in the spatial setting
    coord_names <- coordnames(newdata)
    time_name <- NULL
    
    .plot_spatial_or_ST_common(
        newdata = newdata, column_names = column_names, coord_names = coord_names, 
        time_name = time_name, map_layer = map_layer, 
        palette = palette, plot_over_world = plot_over_world, 
        labels_from_coordnames = labels_from_coordnames, 
        ...
    )
})

#' @rdname plot_spatial_or_ST
#' @export
setMethod("plot_spatial_or_ST", signature(newdata = "SpatialPixelsDataFrame"), 
          function(newdata, column_names,  map_layer = NULL, 
                   subset_time = NULL, palette = "Spectral", 
                   plot_over_world = FALSE, labels_from_coordnames = TRUE,...) {
            
            ## Get the coordinate names, and time_name is NULL in the spatial setting
            coord_names <- coordnames(newdata)
            time_name <- NULL
            
            .plot_spatial_or_ST_common(
              newdata = newdata, column_names = column_names, coord_names = coord_names, 
              time_name = time_name, map_layer = map_layer, 
              palette = palette, plot_over_world = plot_over_world, 
              labels_from_coordnames = labels_from_coordnames, 
              ...
            )
          })


#' @rdname plot_spatial_or_ST
#' @export
setMethod("plot_spatial_or_ST", signature(newdata = "SpatialPolygonsDataFrame"), 
          function(newdata, column_names,  map_layer = NULL, 
                   subset_time = NULL, palette = "Spectral", 
                   plot_over_world = FALSE, labels_from_coordnames = TRUE, ...) {
            
            ## Get the coordinate names, and time_name is NULL in the spatial setting
            coord_names <- coordnames(newdata)
            time_name <- NULL
            
            .plot_spatial_or_ST_common(
              newdata = newdata, column_names = column_names, coord_names = coord_names, 
              time_name = time_name, map_layer = map_layer, 
              palette = palette, plot_over_world = plot_over_world, 
              labels_from_coordnames = labels_from_coordnames, 
              ...
            )
          })

.plot_spatial_or_ST_common <- function(
    newdata, column_names, coord_names, time_name, map_layer, palette, plot_over_world, labels_from_coordnames, ...
) {
  
    ## Inclusion of "-" characters can cause problems; convert to "_"
    column_names <- gsub("-", "_", column_names)
    names(newdata@data) <- gsub("-", "_", names(newdata@data))

    ## Classify the kind of spatial data we have, and convert newdata to a data.frame
    tmp <- .classify_sp_and_convert_to_df(newdata)
    df      <- tmp$df
    sp_type <- tmp$sp_type
    
    ## Remove duplicate columns (can sometimes happen when we convert the Spatial* 
    ## object, e.g., if the coordinates are already present)
    df <- df[, unique(colnames(df))]
    
    # ## Edit the time column so that the facets display t = ...
    # if(!is.null(time_name))
    # df[, time_name] <- factor(
    #    df[, time_name], ordered = TRUE,
    #    labels = paste(time_name, sort(unique(df[, time_name])), sep = " = ")
    # )
    
    if (length(palette) == 1) 
        palette <- rep(palette, length(column_names))
  
    
    ## Plot the requested columns
    plots <- lapply(1:length(column_names), 
                    function(i, x, y, ...) {
                        .plot(df, x[i], coord_names, time_name, sp_type, map_layer, 
                              y[i], plot_over_world, labels_from_coordnames, ...)
                    }, x = column_names, y = palette, ...)
    
    suppressMessages( # suppress message about adding a new coordinate system
        if(plot_over_world) {
            plots <- lapply(plots, function(gg) {
                (gg + .EmptyTheme_theme_only() + coord_map("mollweide")) %>%  
                    draw_world(inc_border = TRUE)
            })
        }
    )
    names(plots) <- column_names
    return(plots)
}

#### This function creates the individual plots.

.plot <- function(df, column_name, coord_names, time_name, sp_type, map_layer, palette, plot_over_world, labels_from_coordnames, ...){
    
    ## Basic plot
    if (!is.null(map_layer)) {
        if ("ggplot" %in% class(map_layer)) {
            ## Do it this way, because we cannot add map_layer to a ggplot object if
            ## map_layer is a gg object itself (e.g., if map_layer is a result from a
            ## call to ggmap(), which is often the case)
            gg <- map_layer 
            gg <- gg %+% df # change the default data to df
        } else {
            ## If map_layer is not a ggplot object, it is probably a layer 
            ## (e.g., the result of a call to geom_polygon)
            gg <- ggplot(df) + map_layer
        }
        
    } else {
        gg <- ggplot(df) 
    }
    ## change/add the default aesthetics, and add some themes for nice plots
    gg <- gg %+% aes_string(x = coord_names[1], y = coord_names[2]) + 
        theme_bw() + coord_fixed()
    
    # ## If NAs are present, tell the user that we hard-code the colour/fill for 
    # ## the pixels with NA values to be transparent 
    # if (any(is.na(df[, column_name]))) {
    #   cat("NA values detected in the data, which will be transparent in the final plot. If you want them to show up in the plot, take the returned plot object and add a colour/fill scale with na.value equal to whatever colour you want.\n ")
    # }
    if (palette == "nasa") {
      colour_fn <- scale_colour_gradientn(colours = nasa_palette, na.value = "transparent")
      fill_fn <- scale_fill_gradientn(colours = nasa_palette, na.value = "transparent")
    } else {
      colour_fn <- scale_colour_distiller(palette = palette, na.value = "transparent")
      fill_fn <- scale_fill_distiller(palette = palette, na.value = "transparent")
    }
    
    ## Plot based on data type
    if (sp_type == "points") { # this is only for observation data
        gg <- gg + geom_point(aes_string(colour = column_name), ...) + colour_fn
    } else {
        if (sp_type == "pixels") {
          ## geom_raster is typically faster, but it doesn't work if we are 
          ## plotting over the plane. For simplicity, just use geom_tile instead.
          if (plot_over_world) {
            gg <- gg + geom_tile(aes_string(fill = column_name), ...)
          } else {
            gg <- gg + geom_raster(aes_string(fill = column_name), ...)
          }
          
        } else if (sp_type == "polygons") {
            gg <- gg + geom_polygon(aes_string(group = "id", fill = column_name), ...) 
        }
        gg <- gg + fill_fn
    }
    
    if (!is.null(time_name))
        gg <- gg + facet_wrap(as.formula(paste("~", time_name)))
    
    if (!labels_from_coordnames)
      gg <- gg + labs(x = expression(s[1]), y = expression(s[2]))

    return(gg)
}


.classify_sp_and_convert_to_df <- function(newdata) {
  
    
    original_names <- names(newdata@data)
    
    if (is(newdata, "SpatialPixelsDataFrame")) {
      df <- cbind(newdata@data, coordinates(newdata))
      sp_type <- "pixels"
    } else if (is(newdata, "SpatialPointsDataFrame")) { ## NB: is(newdata, "SpatialPointsDataFrame") returns TRUE when class(newdata) == "SpatialPixelsDataFrame"
        df <- cbind(newdata@data, coordinates(newdata))
        sp_type <- "points"
    } else if (is(newdata, "SpatialPolygonsDataFrame")) {
        df <- .SpatialPolygonsDataFrame_to_df(newdata)
        sp_type <- "polygons"
    } else if (is(newdata, "STFDF")) {
        if (is(newdata@sp, "SpatialPolygons")) {
            sp_type <- "polygons"
        } else if (is(newdata@sp, "SpatialPixels")) {
            sp_type <- "pixels"
        }
        df <- STFDF_to_df(newdata) 
    } else {
        stop("Class of newdata is not recognised.")
    }

    ## check all original column names are present
    if (!all(original_names %in% names(df)))
      warning("Some of the original names in newdata were not retained when newdata was coerced to a data.frame. 
              This can sometimes happen if some of the column names contain '-', which get converted to '.'")
    
    return(list(df = df, sp_type = sp_type))
}
