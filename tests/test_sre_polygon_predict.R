library(sp)
library(ggplot2)
library(INLA)
library(dplyr)
library(grid)
library(gridExtra)
## Planar SRE

print("Entering test_sre_polygon_predict.R")


test_that("can predict over polygons in plane", {

    ## Get data
    data(meuse)
    meuse$std <- sqrt(0.05066)
    meuse$fs <- 1
    coordinates(meuse) = ~x+y # change into an sp object
    data(meuse.grid)
    gridded(meuse.grid) = ~x + y
    HexPts <- spsample(meuse.grid, type = "hexagonal", cellsize = 200)
    HexPols <- HexPoints2SpatialPolygons(HexPts)
    HexPts_df <- SpatialPointsDataFrame(HexPts,data.frame(fs=rep(1,length(HexPts))))
    HexPols_df <- SpatialPolygonsDataFrame(HexPols,
                                           data.frame(fs=rep(1,length(HexPts)),
                                                      row.names=paste0("ID",rownames(data.frame(HexPts_df)))))
    G <- auto_basis(m = plane(),data=meuse,nres = 2,prune=10,type = "Gaussian")

    test_that("can average basis over polygons in plane", {
        S1 <- eval_basis(G,HexPts_df)
        S2 <- eval_basis(G,HexPols_df)
        plot(as.numeric(S1))
        lines(as.numeric(S2),col='red')
    })

    f <- log(zinc) ~ 1
    S <- SRE(f,meuse,G)
    S <- SRE.fit(S,n_EM = 10,print_lik=T)
    HexPols_df <- SRE.predict(S,HexPols_df,depname = "zinc")
    spplot(HexPols_df,"mu")

})


test_that("can predict over polygons on sphere", {

    df_to_SpatialPolygons <- function(df,keys,coords) {
        df_poly <- plyr::dlply(df,keys,function(d) {
            Sr <- Polygons(list(Polygon(d[coords])),digest::digest(d[keys]))})
        Sr <- SpatialPolygons(df_poly,1:length(df_poly),proj4string=CRS("+proj=longlat"))
    }

    ## Load data
    load(system.file("extdata","AIRS_05_2003.rda", package = "FRK"))
    AIRS_05_2003 <- filter(AIRS_05_2003,day %in% 1) %>%
        mutate(std = co2std, fs = 1,Nobs = 1)
    coordinates(AIRS_05_2003) = ~lon+lat # change into an sp object
    proj4string(AIRS_05_2003)=CRS("+proj=longlat")

    ## Grid prediction
    load(system.file("extdata","isea3h.rda", package = "FRK"))
    ## Need to do res == 5 to start seeing similarity between the point S and the polygon S
    isea3h_res <- filter(isea3h,res == 3) %>%
        arrange(id)
    isea3h_pts <- SpatialPointsDataFrame(filter(isea3h_res,centroid==1)[c("lon","lat")],
                                         data=data.frame(fs=rep(1,sum(isea3h_res$centroid))))

    isea3h_sp <- df_to_SpatialPolygons(df=filter(isea3h_res,centroid==0),keys=c("res","id"),coords=c("lon","lat"))

    isea3h_sp_df <- SpatialPolygonsDataFrame(isea3h_sp,
                                             cbind(filter(isea3h_res,centroid==1) %>% select(-centroid),
                                                   data.frame(fs=rep(1,length(isea3h_sp)),
                                                        row.names=names(isea3h_sp))))

    ## Set up SRE model
    G <- auto_basis(m = sphere(),data=AIRS_05_2003,nres = 2,prune=15,type = "bisquare")
    test_that("can average basis over polygons in sphere", {
        S1 <- eval_basis(G,isea3h_pts)
        S2 <- eval_basis(G,isea3h_sp_df)
    })

    f <- co2avgret ~ lat + 1
    S <- SRE(f,AIRS_05_2003,G)
    S <- SRE.fit(S,n_EM = 3,print_lik=T)
    isea3h_sp_df <- SRE.predict(S,isea3h_sp_df,depname="co2avgret")
#     isea3h_clipped <- left_join(isea3h_res,isea3h_sp_df@data[c("id","mu","var")]) %>%
#                        group_by(res,id) %>%
#                        filter(diff(range(lon)) < 90 & diff(range(lat)) < 90 & centroid==0)

    isea3h_clipped <- left_join(isea3h_res,isea3h_sp_df@data[c("id","mu","var")]) %>%
        filter(centroid ==0) %>%
        clip_polygons_lonlat("id")

    g1 <- (LinePlotTheme() + geom_polygon(data=isea3h_clipped,aes(lon,lat,fill=mu,group=id),colour="light grey") + scale_fill_distiller(palette="Spectral",trans="reverse") + coord_map("mollweide")) %>%
        draw_world(inc_border=FALSE)

    mumin <- min(isea3h_clipped$mu)
    mumax <- max(isea3h_clipped$mu)

    g2 <- (LinePlotTheme() + geom_point(data=data.frame(AIRS_05_2003),aes(lon,lat,colour=pmin(pmax(co2avgret,mumin),mumax)),pch=".") +
               scale_colour_distiller(palette="Spectral",trans="reverse",guide_legend(title="co2")) + coord_map("mollweide")) %>%
        draw_world(inc_border=TRUE)

    #g_all <- arrangeGrob(g1,g2,ncol=1)
})



