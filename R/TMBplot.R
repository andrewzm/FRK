## Plotting functions

# Base it on the names of the supplied data frame 
# i.e. if ("p_mu" %in% colnames(df)) Then plot the mean 
# etc.
# potentially also supply M@response so we can provide discrete scales for bernoulli


#' Plot non-Gaussian data and predictions.
#'
#' @inheritParams .plot_map
#' @inheritParams .plot_data
#' @param zdf A \code{dataframe} containing spatial coordinates (named "x" and "y") and the value of the observations (named "z").
#' @param newdata A \code{dataframe} containing the sspatial coordinates (named "x" and "y"), and the predictions and uncertainty quantification.
#' @return A list of ggplot() objects.
.plot_all <- function(newdata = NULL, zdf = NULL, response = NULL) {
  
  plots <- list()
  
  ## data plot
  if (!is.null(zdf)) plots$data <- .plot_data(zdf, response = response) + labs(colour = 'Z')
  
  ## Prediction and uncertainty plots
  if (!is.null(newdata)) {

    if ("p_Y" %in% names(newdata)) plots$p_Y <- .plot_map(newdata, col = "p_Y") + labs(fill = expression(widehat(p)[Y]["|"][bold(Z)]))
    if ("RMSPE_Y" %in% names(newdata)) plots$RMSPE_Y <- .plot_map(newdata, col = "RMSPE_Y", uncertaintyMap = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[Y]["|"][bold(Z)]))))
    
    if ("p_mu" %in% names(newdata)) plots$p_mu <- .plot_map(newdata, col = "p_mu") + labs(fill = expression(widehat(p)[mu]["|"][bold(Z)]))
    if ("RMSPE_mu" %in% names(newdata)) plots$RMSPE_mu <- .plot_map(newdata, col = "RMSPE_mu", uncertaintyMap = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[mu]["|"][bold(Z)]))))
    
    if ("p_prob" %in% names(newdata)) plots$p_prob <- .plot_map(newdata, col = "p_prob", diverging = TRUE, midpoint = 0.5) + labs(fill = expression(widehat(p)[pi]["|"][bold(Z)]))
    if ("RMSPE_prob" %in% names(newdata)) plots$RMSPE_prob <- .plot_map(newdata, col = "RMSPE_prob", uncertaintyMap = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[pi]["|"][bold(Z)]))))
        
    if ("p_Z_analytic" %in% names(newdata)) plots$p_Z_analytic <- .plot_map(newdata, col = "p_Z_analytic") + labs(fill = expression(widehat(p)[Z]["|"][bold(Z)]~" (analytic)")) 
    if ("p_Z_empirical" %in% names(newdata)) plots$p_Z_empirical <- .plot_map(newdata, col = "p_Z_empirical") + labs(fill = expression(widehat(p)[Z]["|"][bold(Z)]~" (empirical)"))
    if ("RMSPE_Z" %in% names(newdata)) plots$RMSPE_Z <- .plot_map(newdata, col = "RMSPE_Z", uncertaintyMap = TRUE) + labs(fill = expression(sqrt(MSPE(widehat(p)[Z]["|"][bold(Z)]))))
  }

  return(plots)
}


#' Plot of spatial data.
#'
#' @param zdf A \code{dataframe} containing spatial coordinates (named "x" and "y") and the value of the response variable (named "z").
#' @param response A character indicating the assumed response distribution.
#' @param point_size Size of plotted points.
#' @param lim Controls the limits of the colour scale.
#' @return A \code{ggplot} of the data.
#' @seealso \code{\link{.plot_map}}, \code{\link{.plot_all}}
.plot_data <- function(zdf,
                      response,
                      point_size = 1,
                      lim = range(zdf[["z"]])){
  
  ## Make ggplot() object
  p <- ggplot(zdf) +
    geom_point(aes(x, y, colour = z), size = point_size) +
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
#' @param low Low colour (only applicable if \code{diverging == TRUE}).
#' @param mid Mid colour (only applicable if \code{diverging == TRUE}).
#' @param high High colour (only applicable if \code{diverging == TRUE}).
#' @param midpoint Point at which colour \code{mid} is assigned. This is useful for
#' uncertainty ratio maps (in which case a good midpoint is 1) and probability maps
#' (in which case a good midpoint is 0.5).
#' @param uncertaintyMap Logical indicating whether to use uncertainty colour scale.
#' @return A \code{ggplot} of the spatial process.
#' @seealso \code{\link{.plot_data}}, \code{\link{.plot_all}}
.plot_map <- function(df, col, lim = range(df[[col]]), diverging = FALSE,
                     low = "blue", mid = "white", high = "red", midpoint = mean(df[[col]]),
                     uncertaintyMap = FALSE, ...){
  
  ## Basic ggplot() object
  p <- ggplot(df) + geom_tile(aes_string('x','y',fill=col)) +
    theme_bw() + coord_fixed()
  
  ## Colour scale
  if (uncertaintyMap == TRUE & diverging == TRUE) {
    p <- p + scale_fill_gradient2(low = "#01665E", mid = "white", midpoint = midpoint, high = "#8C510A", limits = lim)
  } else if (uncertaintyMap == T) {
    p <- p + scale_fill_distiller(palette = "BrBG", direction = -1, limits = lim)
  } else if (diverging == TRUE) {
    p <- p + scale_fill_gradient2(low = low, mid = mid, midpoint = midpoint, high = high, limits = lim)
  } else{
    p <- p + scale_fill_distiller(palette="Spectral", limits = lim)
  }
  
  return(p)
}






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
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}