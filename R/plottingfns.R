.check_plot_args <- function(args,z=0) {

  if ("max" %in% names(args)) {
    stopifnot(is.numeric(args$max))
  } else {
    args$max <- max(z,na.rm=T)
  }

  if ("min" %in% names(args)) {
    stopifnot(is.numeric(args$min))
  } else {
    args$min <- min(z,na.rm=T)
  }

  if("leg_title" %in% names(args)) {
    stopifnot(is.character(args$leg_title))
  } else {
    args$leg_title <- ""
  }

  if ("plot_dots" %in% names(args)) {
    stopifnot(is.logical(args$plot_dots))
  } else {
    args$plot_dots <- T
  }

  if ("g" %in% names(args)) {
    stopifnot(is.ggplot(args$g))
  } else {
    args$g <- NULL
  }

  if("pt_size" %in% names(args)) {
    stopifnot(is.numeric(args$pt_size))
  } else {
    args$pt_size <- 1
  }

  if("size" %in% names(args)) {
    args$pt_size <- args$size
  }

  if(!("palette" %in% names(args))) {
    args$palette <- "RdYlBu"
  }

  if(!("reverse" %in% names(args))) {
    args$reverse <- FALSE
  }

  return(args)
}

#' @rdname plot_interp
#' @aliases plot_interp,FEBasis,character-method

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}

draw_world <- function(g,inc_border = TRUE) {
    ### Polygon split solution thanks to http://cameron.bracken.bz/finally-an-easy-way-to-fix-the-horizontal-lines-in-ggplot2-maps
    worldmap = map_data("world")
    names(worldmap) <- c("X","Y","PID","POS","region","subregion")
    worldmap = PBSmapping::clipPolys(worldmap, xlim=c(-180,180),ylim=c(-90,90), keepExtra=TRUE)
    if(inc_border) {
        border <- data.frame(X=c(-179.99,-179.99,179.99,179.99),Y=c(-89.99,89.99,89.99,-89.99),PID=1e5,region="border")
        worldmap <- plyr::rbind.fill(worldmap,border)
    }
    g + geom_path(data = worldmap, aes(x=X, y=Y, group=PID), colour="black",size=0.1)
}

setMethod("show_basis",signature(basis = "Basis"),  # GRBF basis with mean offset as last weight
          function(g=ggplot(),basis,s1min = 0, s1max = 1,ds1=0.1,col="red") {
              if(is(manifold(basis),"real_line")) {
                 s <- matrix(seq(s1min,s1max,by=ds1))
                 for (i in 1:basis@n) {
                     S <- basis@fn[[i]](s)
                     df <- data.frame(s=as.numeric(s), y = as.numeric(S))
                     g <- g + geom_line(data=df,aes(x=s,y=y),col=col)
                 }
              } else  if(is(manifold(basis),"plane")) {
                  for (i in 1:basis@n) {
                      g <- g + geom_path(data=circleFun(center=basis@pars[[i]]$loc,
                                                        diameter = basis@pars[[i]]$scale),
                                         aes(x=x,y=y),col=col)
                  }
              } else  if(is(manifold(basis),"sphere")) {
                  for (i in 1:basis@n) {
                      df <-basis@df[,1:2]
                      names(df) <- c("long","lat")
                      g <- g + geom_point(data=df,aes(x=long,y=lat),col=col)
                  }
              }
               return(g)

          })


#' @rdname plot_interp
#' @aliases plot_interp,FEBasis,character-method
setMethod("plot_interp",signature(x = "FEBasis",y="character"),  # GRBF basis with mean offset as last weight
          function(x,y,ds=40,...) {
            args <- list(...)
            Mesh <- x
            df <- getDf(Mesh)
            Mesh['to_plot'] <- df[y]
            args <- .check_plot_args(args,df[y])
            zhi <- args$max
            zlo <- args$min
            leg_title <- args$leg_title
            palette <- args$palette
            reverse <- args$reverse

            mesh_grid <- interp(Mesh["x"],Mesh["y"],Mesh["to_plot"],
                                xo=seq(min(Mesh["x"]),max(Mesh["x"]),length=ds),
                                yo=seq(min(Mesh["y"]),max(Mesh["y"]),length=ds))
            expanded_grid <- cbind(expand.grid(mesh_grid$x,mesh_grid$y),as.vector(mesh_grid$z))
            names(expanded_grid) <- c("x","y","z")
            g <- OverlayPlot(expanded_grid,GGstd = NULL,leg_title=leg_title,zlo=zlo,zhi=zhi,palette=palette,reverse=reverse)
            return(g)
          })

#' @rdname plot_interp
#' @aliases plot_interp,FEBasis-method
setMethod("plot",signature(x = "FEBasis"),  # GRBF basis with mean offset as last weight
          function(x,y,...) {
            args <- list(...)
            args <- .check_plot_args(args)
            plot_dots <- args$plot_dots
            g <- args$g
            g <- ggplotGraph(as.matrix(x@pars$K - diag(diag(x@pars$K))),fixed=x@pars$p,point_size=0.5,g=g,plot_dots=plot_dots)
            return(g)
          })

#' @rdname plot
#' @aliases plot,FEBasis,character-method
setMethod("plot",signature(x = "FEBasis",y="character"),  # GRBF basis with mean offset as last weight
          function(x,y,...) {
            args <- list(...)

            xlocs <- x@pars$p[,1]
            ylocs <- x@pars$p[,2]
            if(!(y %in% names(x@pars$vars))) stop("Did not find label in mesh data frame")
            z <- x@pars$vars[[y]]
            df <- data.frame(x = xlocs, y=ylocs, z=z)
            args <- .check_plot_args(args,z)
            pt_size <- args$pt_size
            df$p1 <- args$max
            df$p2 <- args$min
            leg_title <- args$leg_title

            g <- LinePlotTheme() + geom_point(data=df,aes(x,y,colour=pmax(pmin(z,p1),p2)),size=pt_size) +
              scale_colour_gradient2(low="red",mid="light yellow",high="blue",
                                     guide=guide_legend(title=leg_title,reverse=T))
            return(g)
          })

#' @rdname plot
#' @aliases plot,Obs,character-method
setMethod("plot",signature(x = "Obs",y="character"),  # GRBF basis with mean offset as last weight
          function(x,y,...) {
            args <- list(...)
            if(!(y %in% names(x@df))) stop("Did not find label in mesh data frame")
            z <- x[y]
            args <- .check_plot_args(args,z)
            lab <- args$leg_title
            pt_size <- args$pt_size
            p1 <- args$max
            p2 <- args$min

            xlocs <- x["x"]
            ylocs <- x["y"]

            df <- data.frame(x = xlocs, y=ylocs, z=z)
            df$p1 <- p1
            df$p2 <- p2

            g <- LinePlotTheme() + geom_point(data=df,aes(x,y,colour=pmax(pmin(z,p1),p2)),size=pt_size) +
              #scale_colour_gradient2(low=muted("red"),mid="light yellow",high=muted("blue"),guide=guide_legend(lab))
              scale_colour_distiller(palette="RdYlBu",guide=guide_legend(lab,reverse=T))
            return(g)
          })

#' @rdname plot
#' @aliases plot,Obs_poly,character-method
setMethod("plot",signature(x = "Obs_poly",y="character"),  # GRBF basis with mean offset as last weight
          function(x,y,...) {
            args <- list(...)
            if(!(y %in% names(x@df2))) stop("Did not find label in mesh data frame")
            z <- x[y]
            args <- .check_plot_args(args,z)
            lab <- args$leg_title
            pt_size <- args$pt_size
            p1 <- args$max
            p2 <- args$min

            df <- x@df2
            df$p1 <- p1
            df$p2 <- p2
            df$to_plot <- df[y][,1]

            g <- LinePlotTheme() + geom_polygon(data=df,aes(x=x,y=y,fill=pmax(pmin(to_plot,p1),p2),group=id)) +
              scale_fill_gradient2(low = muted("red"),mid="light yellow",high=muted("blue"),
                                   guide=guide_legend(lab))
            return(g)
          })


## Draw a graph of the precision matrix
DrawGraph <- function(Q,fixed=NA)   {
  diag(Q) <- 0
  Q[which(Q!=0)] = 1 		# make symbolic
  mygraph <- graph.adjacency(Q,mode="undirected")
  if(is.na(fixed)){
    layout =  layout.fruchterman.reingold(mygraph)
  } else {
    layout =  fixed
  }


  if (dim(Q)[1] < 200) {
    plot(mygraph,layout=layout,vertex.size=4,vertex.label.dist=0.5,vertex.color="red",edge.arrow.size=0.5)
  } else {
    plot(mygraph,layout=layout,vertex.size=4,vertex.label.dist=NA,vertex.color="red",edge.arrow.size=0.5)
  }
}



plot_scale_bar <- function(g_orig,x=0,y=0,l=10)  {
  df1 <- data.frame(x = c(x, x = x + l/2,  x = x + l/2,  x = x),
                    y = c(y + l/20,y + l/20, y - l/20, y - l/20))
  df2 <- data.frame(x = c(x, x = x - l/2,  x = x - l/2,  x = x),
                    y = c(y + l/20,y + l/20, y - l/20, y - l/20))
  dftext <- data.frame(x = x, y = y - l/3, l= l)
  g <- g_orig + geom_polygon(data=df1,aes(x,y),colour="black",fill="white") +
    geom_polygon(data=df2,aes(x,y),colour="black",fill="black") +
    geom_text(data=dftext,aes(x,y,label = paste(l,'km')))
  return(g)

}

ggplotGraph <- function(Q,fixed=NA,point_size=0.5,g = NULL,plot_dots=T) {

  net <- network(Q,directed = F,density=0.3)
  m <- as.matrix.network.adjacency(net) # get sociomatrix
  plotcord <- data.frame(X1=fixed[,1],X2=fixed[,2])
  edglist <- as.matrix.network.edgelist(net)
  edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
  colnames(edges) <-  c("X1","Y1","X2","Y2")
  edges$midX  <- (edges$X1 + edges$X2) / 2
  edges$midY  <- (edges$Y1 + edges$Y2) / 2
  if (class(g)[1] == "NULL") {
    pnet <- ggplot()  +
      geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),
                   data=edges, size = 0.3, colour="dark grey") +
      scale_colour_brewer(palette="Set1") +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
      # discard default grid + titles in ggplot2
      theme(panel.background = element_blank()) + theme(legend.position="none")+
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      theme( legend.background = element_rect(colour = NA)) +
      theme(panel.background = element_rect(fill = "white", colour = NA)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  } else {
    pnet <- g  +
      geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),
                   data=edges, size = 0.3, colour="dark grey") +
      scale_colour_brewer(palette="Set1")
    #scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL)
    # discard default grid + titles in ggplot2
    #theme(panel.background = element_blank(), legend.position="none",
    #      legend.background = element_rect(colour = NA), panel.background = element_rect(fill = "white", colour = NA),
    #      panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  }
  if(plot_dots) pnet <- pnet + geom_point(aes(X1, X2),colour="red", data=plotcord,size=point_size)
  return(pnet)

}

## Plot a spatial field from vertices
PlotSpatField <- function(p,x,my_title = '',ylab="",zlim=c(0,0),palette=brewer.pal(11,"RdBu"),res=200)  {
  surf <- interp(x=p[,1],y=p[,2],z=x,xo=seq(min(p[,1]), max(p[,1]), length = res),yo=seq(min(p[,2]), max(p[,2]), length = res))
  if ((zlim[1] == 0) && (zlim[2] == 0)) zlim <- range(surf$z,na.rm=T)
  zlen <- zlim[2] - zlim[1]

  surf$z[which(surf$z < zlim[1])] = zlim[1]
  surf$z[which(surf$z > zlim[2])] = zlim[2]


  image(x=surf$x, y=surf$y, z=surf$z,col=palette, axes=F, zlim=zlim, xlab="", ylab=ylab)
  image.plot( zlim=zlim, col=palette, legend.only=TRUE, horizontal =TRUE)
  grid()
  box()
  # filled.contour(surf$x, surf$y, surf$z,color.palette=terrain.colors,xlab='x',ylab='y')
  title(my_title)
  return(surf)
}

## Plot from Voronoi tesselations
DrawfromVoronoi <- function(Voronoi,col,p) {
  n = length(col)
  for (i in 1:n)  {
    edges <- rbind(as.matrix(subset(Voronoi$dirsgs,ind1 ==i,select=c(x1,y1,bp1))),
                   as.matrix(subset(Voronoi$dirsgs,ind2 ==i,select=c(x1,y1,bp1))),
                   as.matrix(subset(Voronoi$dirsgs,ind1 ==i,select=c(x2,y2,bp2))),
                   as.matrix(subset(Voronoi$dirsgs,ind2 ==i,select=c(x2,y2,bp2))))
    X <- unique(edges)
    if(sum(X[,3])>0) {
      X <- rbind(X,c(p[i,],0))
    }
    edges <- X[,(1:2)]
    edges <- edges[chull(edges), ]
    polygon(unname(edges[,1]),unname(edges[,2]),col=col[i],border='NA')
  }
}

#' @title Line-plot theme
#'
#' @description Formats a ggplot object for neat plotting.
#'
#' @return Object of class \code{ggplot}
#' @keywords ggplot
#' @export
#' @examples
#' \dontrun{
#' X <- data.frame(x=runif(100),y = runif(100), z = runif(100))
#' LinePlotTheme() + geom_point(data=X,aes(x,y,colour=z))
#'}
LinePlotTheme <- function() {
  g <- ggplot() + theme(panel.background = element_rect(fill='white', colour='black'),text = element_text(size=20),
                        panel.grid.major =  element_line(colour = "light gray", size = 0.05),
                        panel.border  = element_rect(fill=NA, colour='black'),
                        plot.margin=unit(c(5,5,5,0),"mm"))
  return (g)
}

#' @title Empty-plot theme
#'
#' @description Formats a ggplot object for plotting with no annotations/grids.
#'
#' @return Object of class \code{ggplot}
#' @keywords ggplot
#' @export
#' @examples
#' \dontrun{
#' X <- data.frame(x=runif(100),y = runif(100), z = runif(100))
#' EmptyTheme() + geom_point(data=X,aes(x,y,colour=z))
#'}
EmptyTheme <- function() {
  g <- ggplot() +  theme(panel.background = element_rect(fill='white', colour='white'),panel.grid=element_blank(),axis.ticks=element_blank(),
                         panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
  return (g)
}
## this temperature field suggested in Kaye et al (2011), Fig 3 in the
## response to my comment
makeKaye1 <- function(mbr = c(-15, 0, 10, 25), sdbr = c(0.3, 0.6, 1),rtb=F) {

  stopifnot(length(mbr) == 4, length(sdbr) == 3)

  if (rtb == F) {
    sat <- cbind(red = c(0, 100, 240, 230, 170),
                 green = c(0, 180, 210, 120, 0),
                 blue = c(170, 255, 0, 0, 0))
  } else {
    sat <- cbind(red = c(170,230,240,100,0),
                 green = c(0,120,210,180,0),
                 blue = c(0,0,0,255,170))
  }


  pp <- c(0, 0.33, 0.66, 1)
  coltab <- sapply(pp, function(x) {
    y <- (1 - x) * sat + x * 255
    rgb(y[, "red"], y[, "green"], y[, "blue"], alpha = 255, maxColorValue = 255)
  })

  rownames(coltab) <- c(paste("<", mbr[1]),
                        paste(mbr[-length(mbr)], mbr[-1], sep = " to "),
                        paste(">=", mbr[length(mbr)]))
  mbr <- c(-Inf, mbr, Inf) # for cut

  sdbr <- c(0, sdbr)
  colnames(coltab) <- c(paste(sdbr[-length(sdbr)], sdbr[-1], sep = " to "),
                        paste(">=", sdbr[length(sdbr)]))
  sdbr <- c(sdbr, Inf) # for cut

  function(m, sd = 0) {
    m <- cut(m, breaks = mbr, right = FALSE, labels = FALSE) # in 1:5 (small - large)
    sd <- cut(sd, breaks = sdbr, right = FALSE, labels = FALSE)  # in 1:4 (unc - certain)
    coltab[cbind(m, sd)]
  }
}
mymap <- function(grid, col, Dx,
                  xlim = c(-180, 180), ylim = c(-85, 85),
                  add = FALSE,
                  coast = c("under", "over", "none"), colcoast = "darkgrey",
                  boxtype = c("ticks", "none", "box")) {

  ## simple checks

  stopifnot(is.data.frame(grid), c("lon", "lat") %in% names(grid),
            nrow(grid) == length(col))
  coast <- match.arg(coast)
  boxtype <- match.arg(boxtype)

  ## fill in gridcell widths

  if (missing(Dx)) {
    Dx <- sapply(grid[, c("lon", "lat")], function(x)
      min(diff(unique(sort(x)))))
    #     cat(sprintf("** \'Dx\' inferred from grid, dlon = %.2f, dlat = %.2f\n",
    #                 Dx["lon"], Dx["lat"]))
  }

  ## initialise plot and plot gridcells

  #if (!add)
  #  map(xlim = xlim, ylim = ylim, type = "n", mar = par("mar"))
  #if (coast == "under")
  #  map(interior = FALSE, add = TRUE, col = colcoast, mar = par("mar"))

  rect(grid$lon, grid$lat, grid$lon + Dx["lon"], grid$lat + Dx["lat"],
       col = col, border = col, lwd = 0, ljoin = 1)

  #if (coast == "over")
  #  map(interior = FALSE, add = TRUE, col = colcoast, mar = par("mar"))

  if (boxtype == "box") box()
  else if (boxtype == "ticks") {
    axis(1, at = pretty(xlim, 10), tcl = 0.25, cex.axis = 0.7,
         mgp = c(2, 0.0, 0))
    axis(2, at = pretty(ylim, 7), tcl = 0.25, cex.axis = 0.7,
         mgp = c(2, 0.25, 0), las = 1)
    box()
  }

  invisible(NULL)
}
DrawKaylePalette <- function(coltab)
{
  ## create legend as separate figure
  do.call(getOption("device"), list(width = 4, height = 3.5))
  coltab <- get("coltab", envir = environment(Kaye))
  par(mar = c(4, 5.5, 1, 1), mgp = c(1.25, 0.25, 0), las = 1, cex = 0.8)

  plot.new()
  plot.window(xlim = c(0, 4), ylim = c(0, 5))

  xleft <- rep(0:3, 5); ybottom = rep(0:4, rep(4, 5))
  rect(xleft, ybottom, xleft + 1, ybottom + 1, col = t(coltab), border = "darkgrey")
  axis(1, at = 0:3 + 0.5, colnames(coltab), tick = FALSE)
  #title(xlab = "Standard deviation (m/a)", line = 2)
  axis(2, at = 0:4 + 0.5, rownames(coltab), tick = FALSE)
  #title(ylab = "Height change (m/a)", line = 3.75)

}
## Draw a border using countours from a mask of 1/0s
DrawBorderFromMask <- function(x,y,z)    {
  contour(unique(x),unique(y),z,axis=T,levels=seq(from=-3, to=3, by=1),drawlabels = F,add=T,lwd=1)
}
## Plot gridded data from N x 1 vectors
Plot_gridded_data <- function(x,y,z,palette=rainbow(128),zlim=c(-1,1),xlab="",ylab="")    {
  x <- unique(x)
  y <- unique(y)
  z <- reshape(matrix(z),length(y),length(x))
  image(x=x, y=y, z=t(z),col=palette, axes=F, zlim=zlim, xlab=xlab, ylab=ylab)
  image.plot( zlim=zlim, col=palette, legend.only=TRUE, horizontal =TRUE)
  grid()
  box()
}
Plotsigmafromkappa <- function(sigma,a,b) {
  plot(sigma,2*sigma^(-3)*dgamma(1/(sigma^2),shape=a,rate=b),ylab="")
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
get_ellipse <- function(data, fill){
  edata <- as.matrix(data)
  ehull <- ellipsoidhull(edata)
  phull <- as.data.frame(predict(ehull))
  data.frame(
    x=phull$V1,
    y=phull$y,
    fill=rep(fill, nrow(phull))
  )
}
spy <-function(X) {
  Y <- which(t(X)!=0,arr.ind=TRUE)
  plot(Y[,1],-Y[,2],pch=".")
}

#' @title Overla Plot
#' @description Still to be documented
#' @export
OverlayPlot <- function(GG,GGstd = NULL,leg_title="",do_mesh=0,zlo = -0.6,zhi = 0.6,alphalo = 0.2, alphahi = 1,palette="RdYlBu",reverse=F) {
  GG$zlo = zlo
  GG$zhi = zhi
  if (!is.null(GGstd)) {
    GGstd$alphalo = alphalo
    GGstd$alphahi = alphahi
  }

  g <- LinePlotTheme() + geom_tile(data=subset(GG,!(is.na(z))),aes(x=x,y=y,fill=pmin(pmax(z,zlo),zhi))) +
    #scale_fill_gradient2(low=(muted("red")),mid="light yellow",high=(muted("blue")),limits=c(zlo,zhi),
    #                     guide=guide_colourbar(title=leg_title)) +
    scale_fill_distiller(palette=palette,reverse=reverse,limits=c(zlo,zhi),
                         guide=guide_colourbar(title=leg_title,position=c(-1500,500)))

  if (!is.null(GGstd)) {
    g <- g + geom_tile(data=GGstd,aes(x=x,y=y, alpha=-pmax(pmin(z,alphahi),alphalo)),fill="white") +scale_alpha(guide = 'none')
  }


  if (do_mesh) g <- ggplotGraph(as.matrix(K_ICE - diag(diag(K_ICE))),fixed=p_ICE,point_size=0.01,g=g,plot_dots=F)
  g <- g + theme(legend.position="right") +
    theme(legend.key.width=unit(4,"lines")) + theme(legend.key.height=unit(4,"lines")) +
    theme(text = element_text(size=40))
  return(g)
}

# colourful palette for detailed plots
RdGrPu <- function(int=3) {

  c(brewer_pal("seq", "Reds")(int), brewer_pal("seq", "Oranges")(int),
    brewer_pal("seq", "Greens")(int),brewer_pal("seq", "Greys")(int),
    brewer_pal("seq", "Purples")(int),brewer_pal("seq", "Blues")(int))

}

medal_plot <- function(obs,Qx,inc.matrix,g=ggplot(),print_middle=T,alpha=1,scale=0.004,Sigma_part=NULL,xsamp=NULL,clamp_below=F) {
  stopifnot(!(is.null(Sigma_part)) | !(is.null(xsamp)))

  C_full <- inc.matrix
  cholQx <- cholPermute(Qx)
  Qobs <- diag(1/(obs$std^2))

  PP <- t(cholsolve(Qx,t(C_full),perm=T,cholQp = cholQx$Qpermchol, P = cholQx$P))  # A backsolve is a transposed forward solve ...
  Prior_var <- rowSums(PP * C_full)

  PP <- t(cholsolve(Qx,t(C_full),perm=T,cholQp = cholQx$Qpermchol, P = cholQx$P))  # A backsolve is a transposed forward solve ...
  Prior_var_diag <- diag(diag((PP %*% t(C_full))))
  Obs_var <- 1/diag(Qobs)
  Var_tot2 <- diag(solve(solve(Prior_var_diag) + solve(diag(Obs_var))))


  if(!(is.null(Sigma_part))) {

    PP <- C_full %*% Sigma_part
    Var_tot <- rowSums(PP * C_full)


    #dev.new(); plot(Var_tot/Obs_var)
    #dev.new(); plot(Var_tot/Prior_var)
  }   else {
    Var_tot <- rep(0,nrow(C_full))
    for (i in 1:nrow(C_full)) {
    #Sigma <- (x_samp[,((burnin+1):m)] - x_des$mean) %*% t(x_samp[,((burnin+1):m)] - x_des$mean)/(m - burnin)
    ind <- which(!(C_full[i,] == 0))
    Cov <- (x_samp[ind,] - apply(x_samp[ind,],1,mean)) %*% t(x_samp[ind,] - apply(x_samp[ind,],1,mean))/(ncol(xsamp)-1)
    Var_tot[i] <- C_full[i,ind] %*% Cov %*% (C_full[i,ind])
    }
  }
  obs$PostVar <- Var_tot
  obs$PostVar2 <- Var_tot2
  obs$ObsVar <- Obs_var
  obs$PriorVar <- Prior_var
  obs$Outer_disk_val <- pmin(Obs_var,Prior_var)
  #obs$Outer_disk_val <- Prior_var
  obs$Outer_disk_col <- apply(cbind(Obs_var,Prior_var),1,function(x) {
    minbnd <- which.min(x)
    if(minbnd == 1) {
      return("#3A3A98FF")
    } else {
      return("#832424FF")
    }})


  size_outer <-  scale*(sqrt(obs$Outer_disk_val))
  size_inner <- scale*(sqrt(obs$PostVar))
  size_middle <- scale*(sqrt(obs$PostVar2))
  size_middle <- pmin(size_middle,0.9*size_outer) # this is to show the rim colour
  size_inner <- pmin(size_inner,0.9*size_middle)

  ## If medals are too small scale them up
  if(clamp_below) {
    min_size <- min(diff(range(obs$x)),diff(range(obs$y)))/900
    ind_too_small <- which(size_outer < min_size)   # don't make medals too small
    scales <- min_size/size_outer
    size_outer[ind_too_small] <- min_size   # don't make medals too small
    size_inner[ind_too_small] <- (size_inner*scales)[ind_too_small]
    size_middle[ind_too_small] <- (size_middle*scales)[ind_too_small]
  }

  mobs <- nrow(obs)
  for (i in 1:mobs) {
    g <- g + geom_polygon(data=circleFun(c(obs$x[i],obs$y[i]),15*size_outer[i],npoints = 100),aes(x,y),fill=obs$Outer_disk_col[i],alpha=alpha)
  }
  if(print_middle)
    for (i in 1:mobs) {
      g <- g + geom_polygon(data=circleFun(c(obs$x[i],obs$y[i]),15*size_middle[i],npoints = 100),aes(x,y),fill="white",alpha=alpha)
    }
  for (i in 1:mobs) {
    g <- g + geom_polygon(data=circleFun(c(obs$x[i],obs$y[i]),15*size_inner[i],npoints = 100),aes(x,y),fill="#FFD700") #FF8800
  }
  return(g)
}

