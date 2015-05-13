library(sp)
library(ggplot2)
library(INLA)
library(dplyr)
library(grid)
library(gridExtra)

print("Entering test_sre_point_predict.R")

## Planar SRE
test_that("can do point predictions using planar SRE", {

    ## Get data
    data(meuse)
    meuse$std <- sqrt(0.05066)/5
    meuse$fs <- 1
    coordinates(meuse) = ~x+y # change into an sp object

    ## Construct SRE
    f <- log(zinc) ~ 1
    G <- auto_basis(m = plane(),data=meuse,nres = 2,prune=10,type = "Gaussian")
    S <- SRE(f,meuse,G)
    S <- SRE.fit(S,n_EM = 10,print_lik=T)

    ## Prediction grid and prediction.
    xrange <- range(meuse$x)
    yrange <- range(meuse$y)
    xgrid <- seq(xrange[1],xrange[2],by=50)
    ygrid <- seq(yrange[1],yrange[2],by=50)
    xy <- expand.grid(x=xgrid,y=ygrid)
    border <- INLA::inla.nonconvex.hull(coordinates(meuse),convex = -0.05)$loc
    pred_locs <- subset(SDMTools::pnt.in.poly(xy,border),pip==1)
    pred_locs$fs <- 1
    pred_locs$pip <- NULL
    coordinates(pred_locs) <- ~x + y
    pred_locs <- SRE.predict(S,pred_locs,depname = "zinc")

    ## Plot
    g <- LinePlotTheme() + geom_tile(data=data.frame(pred_locs),aes(x,y,fill=mu)) + geom_point(data=data.frame(meuse),aes(x,y,col=log(zinc)),size=3)  + scale_colour_distiller(palette = "Spectral") + scale_fill_distiller(palette="Spectral") + coord_fixed()
    show_basis(g,S@basis)
})

## Spherical SRE
test_that("can do point predictions using spherical SRE", {

    ## Load data
    load(system.file("extdata","AIRS_05_2003.rda", package = "FRK"))
    AIRS_05_2003 <- filter(AIRS_05_2003,day %in% 1) %>%
                    mutate(std = co2std, fs = 1)
    coordinates(AIRS_05_2003) = ~lon+lat # change into an sp object

    ## Set up SRE model
    f <- co2avgret ~ lat + 1
    G <- auto_basis(m = sphere(),data=AIRS_05_2003,nres = 2,prune=15,type = "bisquare")
    S <- SRE(f,AIRS_05_2003,G)

    ## Fit model
    S <- SRE.fit(S,n_EM = 3,print_lik=T)

    ## Prediction locations
    load(system.file("extdata","isea3h.rda", package = "FRK"))
    isea3h_c <- filter(isea3h,res == 6 & centroid==1) %>%
                mutate(fs = 1)
    coordinates(isea3h_c) <- ~lon + lat
    isea3h_c <- SRE.predict(S,isea3h_c,depname="co2avgret")

    ## Convert to data frame and clip for plotting
    isea3h_clipped <- isea3h %>%
              filter(res==6) %>%
              left_join(data.frame(isea3h_c)[c("id","res","mu","var")]) %>%
                filter(centroid ==0) %>%
                clip_polygons_lonlat("id")


    ### Implement a clipping procedure... which makes two polygons out of one at the boundary. Maybe
    ### see the code of PBSmapping::clipPolys

    g1 <- (LinePlotTheme() + geom_polygon(data=isea3h_clipped,aes(lon,lat,fill=mu,group=id),colour="light grey") + scale_fill_distiller(palette="Spectral",trans="reverse") + coord_map("mollweide")) %>%
        draw_world(inc_border=FALSE)

    mumin <- min(isea3h_clipped$mu)
    mumax <- max(isea3h_clipped$mu)

    g2 <- (LinePlotTheme() + geom_point(data=data.frame(AIRS_05_2003),aes(lon,lat,colour=pmin(pmax(co2avgret,mumin),mumax)),pch=".") +
               scale_colour_distiller(palette="Spectral",trans="reverse",guide_legend(title="co2")) + coord_map("mollweide")) %>%
        draw_world(inc_border=TRUE)

    g3 <- (LinePlotTheme() + geom_polygon(data=isea3h_clipped,aes(lon,lat,fill=pmin(sqrt(var),2.74),group=id),colour="light grey") + scale_fill_distiller(palette="Spectral",trans="reverse") + coord_map("mollweide")) %>%
        draw_world(inc_border=FALSE)

    #g_all <- arrangeGrob(g1,g2,ncol=1)

})
