#' @title Draw a map of the world with country boundaries.
#' @description Generates a ggplot2 map of the world
#' @param g initial ggplot object
#' @param inc_border flag indicating whether a map border should be drawn or not; see details.
#' @details This function uses \code{ggplot2::map_data} in order to create a world map. Since, by default, this creates lines crossing the world at the (-180,180) longitude boundary, function \code{clipPolys} taken from the package \code{PBSmapping} is used to split the polygons at this boundary into two. If \code{inc_border} is TRUE, then a border is drawn around the lon/lat space; this option is most useful for projections which do not yield rectangular plots (e.g., the sinusoidal projection).
#' @export
#' @examples
#' library(ggplot2)
#' draw_world(g = ggplot())
draw_world <- function(g = ggplot(),inc_border = TRUE) {
    ### Polygon split solution thanks to http://cameron.bracken.bz/finally-an-easy-way-to-fix-the-horizontal-lines-in-ggplot2-maps
    if(!(is(g, "ggplot"))) stop("g has to be of class ggplot")
    if(!(is.logical(inc_border))) stop("inc_border needs to be TRUE or FALSE")
    X <- Y <- PID <- NULL
    data(worldmap, envir=environment(), package = "FRK")
    names(worldmap) <- c("X","Y","PID","POS","region","subregion")
    worldmap = clipPolys(worldmap, xlim=c(-180,180),ylim=c(-90,90), keepExtra=TRUE)
    if(inc_border) {
        border <- data.frame(X=c(-179.99,-179.99,179.99,179.99),Y=c(-89.99,89.99,89.99,-89.99),PID=1e5,region="border")
        worldmap <- plyr::rbind.fill(worldmap,border)
    }
    g + geom_path(data = worldmap, aes(x=X, y=Y, group=PID), colour="black",size=0.1)
}

#' @rdname show_basis
#' @aliases show_basis,Basis-method
setMethod("show_basis",signature(basis = "Basis"),  # GRBF basis with mean offset as last weight
          function(basis,g=ggplot()) {

              message("Note: show_basis assumes spherical distance functions when plotting")

              y <- res <- x <- lon <- lat <- NULL # Suppress bindings warning

              if(is(manifold(basis),"real_line")) {
                 s1min <- min(basis@df$loc1)  - max(basis@df$scale)*3
                 s1max <- max(basis@df$loc1)  + max(basis@df$scale)*3
                 s <- matrix(seq(s1min,s1max,length=1000))
                 for (i in 1:basis@n) {
                     S <- basis@fn[[i]](s)
                     df <- data.frame(s=as.numeric(s), y = as.numeric(S),res=basis@df$res[i])
                     g <- g + geom_line(data=df,aes(x=s,y=y,col=as.factor(res))) + labs(colour="res")
                 }
              } else  if(is(manifold(basis),"plane")) {
                  for (i in 1:basis@n) {
                      g <- g + geom_path(data=cbind(circleFun(center=as.numeric(basis@df[i,1:2]),
                                                        diameter = basis@df$scale[i]),
                                                    res=as.factor(basis@df$res[i])),
                                         aes(x=x,y=y,linetype=res))
                  }
              } else  if(is(manifold(basis),"sphere")) {
                  df <-basis@df
                  df <- df[rev(rownames(df)),]
                  names(df)[1:2] <- c("lon","lat")
                  g <- g + geom_point(data=df,aes(x=lon,y=lat,size=res),shape=1) +
                           scale_size_continuous(trans="reverse",breaks =1:10)

              } else  if(is(manifold(basis),"STplane")) {
                  df <-basis@df
                  df <- df[rev(rownames(df)),]
                  names(df)[1:2] <- c("x","y")
                  g <- g + geom_point(data=df,aes(x=x,y=y,size=res),shape=1) +
                      scale_size_continuous(trans="reverse",breaks =1:10) +
                      facet_wrap(~loc3)

              } else  if(is(manifold(basis),"STsphere")) {
                  df <-basis@df
                  df <- df[rev(rownames(df)),]
                  names(df)[1:2] <- c("lon","lat")
                  g <- g + geom_point(data=df,aes(x=lon,y=lat,size=res),shape=1) +
                      scale_size_continuous(trans="reverse",breaks =1:10) +
                      facet_wrap(~loc3)

              }
               return(g)

          })


#' @name plotting-themes
#' @aliases LinePlotTheme
#' @aliases EmptyTheme
#' @title Plotting themes
#' @description Formats a ggplot object for neat plotting.
#' @return Object of class \code{ggplot}
#' @keywords ggplot
#' @export
#' @details \code{LinePlotTheme()} creates \code{ggplot} object with a white background, a relatively large font and grid lines. \code{EmptyTheme()} on the other hand creates a \code{ggplot} object with no axes or legends.
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
                        #plot.margin=unit(c(5,5,5,0),"mm"))
  return (g)
}

#' @rdname plotting-themes
#' @export
EmptyTheme <- function() {
  g <- ggplot() +  theme(panel.background = element_rect(fill='white', colour='white'),panel.grid=element_blank(),axis.ticks=element_blank(),
                         panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
  return (g)
}


clip_polygons_lonlat <- function(d,key) {
    lon <- lat <- NULL
    plyr::ddply(d,key,function(df) {
        if(diff(range(df$lon)) > 90) {
            Y1 <- filter(df,lon >= 0) %>% mutate(id= df[key][1,]*1e6)
            Y1$lon[which(Y1$lon %in% sort(Y1$lon,decreasing=T)[1:2])] <- 179.99
            Y1 <- rbind(Y1,Y1[1,])
            Y2 <- filter(df,lon <= 0) %>% mutate(id= df[key][1,]*1e6+1)
            Y2$lon[which(Y2$lon %in% sort(Y2$lon,decreasing=F)[1:2])] <- -179.99
            Y2 <- rbind(Y2,Y2[1,])
            rbind(Y1,Y2)
        } else {
            df
        }})
}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}
