## Plotting functions


# Potentially supply M@response so we can provide discrete scales for bernoulli


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





## I think it is better to use a dedicated plotting package for this purpose, 
## such as cowplot or ggpubr.
## Aligning ggplots
.align_plots <- function(...) {

  if (!requireNamespace(c("grid", "gtable"), quietly = TRUE)) {
    stop("Packages \"grid\" and \"gtable\" needed for this function to work. Please install them.",
         call. = FALSE)
  }

  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]

  plots.grobs <- lapply(list(...), ggplotGrob)

  max.widths <- do.call(grid::unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })

  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (gtable::is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable::gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })

  plots.grobs.eq.widths.aligned
}