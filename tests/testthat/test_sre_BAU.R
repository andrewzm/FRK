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
    meuse$fs <- 1
    coordinates(meuse) = ~x+y # change into an sp object
    data(meuse.grid)
    gridded(meuse.grid) = ~x + y
    HexPts <- spsample(meuse.grid, type = "hexagonal", cellsize = 200)
    HexPols <- HexPoints2SpatialPolygons(HexPts)
    HexPts_df <- SpatialPointsDataFrame(HexPts,data.frame(fs=rep(1,length(HexPts))))
#     HexPols_df <- SpatialPolygonsDataFrame(HexPols,
#                                            data.frame(fs=rep(1,length(HexPts)),
#                                                       row.names=paste0("ID",rownames(data.frame(HexPts_df)))))
    HexPols_df <- aggregate(meuse,HexPols)
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

    map_data_to_BAUs <- function(data_sp_pts,sp_pols,av_var,variogram.formula=NULL) {
        data_sp_pts$Nobs <- 1
        Data_in_BAU <- over(sp_pols,data_sp_pts[c(av_var,"Nobs")],fn=sum)
        sp_pols@data[av_var] <- Data_in_BAU[av_var]/Data_in_BAU$Nobs
        sp_pols@data["Nobs"] <- Data_in_BAU$Nobs

        new_sp_pts <- SpatialPointsDataFrame(coords=sp_pols[coordnames(data_sp_pts)]@data,
                                             data=sp_pols@data,
                                             proj4string = CRS(proj4string(data_sp_pts)))
        new_sp_pts$std <- sqrt(1 / new_sp_pts$Nobs)
        new_sp_pts <- subset(new_sp_pts,!is.na(Nobs))

        if(is(variogram.formula,"formula")) {
            v <- gstat::variogram(object=variogram.formula,#co2avgret~lat+1,
                                        data=new_sp_pts,cressie=T)
            warning("Not accounting for multiple data in the same grid box during variogram estimation. Need to see how to do this with gstat")
            #AIRS.vgm = gstat::variogram(co2avgret~lat+1, AIRS_05_2003,cressie=T)
            vgm.fit = gstat::fit.variogram(v, model = gstat::vgm(1, "Sph", 2000, 1))
            plot(v,vgm.fit)
            print(paste0("sigma2e estimate = ",vgm.fit$psill[1]))
            new_sp_pts$std <- sqrt(vgm.fit$psill[1] / new_sp_pts$Nobs)
        }
        new_sp_pts
    }



    ## Load data
    load(system.file("extdata","AIRS_05_2003.rda", package = "FRK"))
    AIRS_05_2003 <- filter(AIRS_05_2003,day %in% 1) %>%
        select(lon,lat,co2avgret)
    coordinates(AIRS_05_2003) = ~lon+lat # change into an sp object
    proj4string(AIRS_05_2003)=CRS("+proj=longlat")

    ## Prediction (BAU) grid
    load(system.file("extdata","isea3h.rda", package = "FRK"))
    isea3h_res <- filter(isea3h,res == 6) %>%
        arrange(id) %>%
        group_by(id) %>%
        filter(diff(range(lon)) < 90) %>% data.frame()
    isea3h_res_centroids <- filter(isea3h_res,centroid==1)
    isea3h_res_centroids$fs <- 1

    isea3h_sp_pts <- SpatialPointsDataFrame(isea3h_res_centroids[c("lon","lat")],
                                         data=data.frame(fs=rep(1,sum(isea3h_res$centroid)),
                                                         id = isea3h_res_centroids$id),
                                         proj4string=CRS("+proj=longlat"))

    isea3h_sp_pol <- df_to_SpatialPolygons(df=filter(isea3h_res,centroid==0),
                                           keys=c("id"),
                                           coords=c("lon","lat"))

#     isea3h_sp_poldf <- SpatialPolygonsDataFrame(isea3h_sp_pol,
#                                             cbind(data.frame(row.names=names(isea3h_sp_pol)),
#                                                   coordinates(isea3h_sp_pts),
#                                                   id=isea3h_sp_pts@data$id))

    isea3h_sp_poldf <- SpatialPolygonsDataFrame(isea3h_sp_pol,
                                                cbind(data.frame(row.names=names(isea3h_sp_pol)),
                                                      isea3h_res_centroids[c("id","lon","lat","fs")]))

    AIRS_remapped <- map_data_to_BAUs(AIRS_05_2003,
                                      isea3h_sp_poldf,
                                      av_var = "co2avgret",
                                      variogram.formula = co2avgret~lat+1)
    AIRS_remapped$fs <- 1

    ## Set up SRE model
    G <- auto_basis(m = sphere(),data=AIRS_remapped,nres = 3,prune=15,type = "bisquare")
    f <- co2avgret ~ lat + 1
    S <- SRE(f,AIRS_remapped,G)
    S <- SRE.fit(S,n_EM = 3,print_lik=T)
    isea3h_sp_pts <- SRE.predict(S,isea3h_sp_pts,depname="co2avgret")

    isea3h_clipped <- left_join(isea3h_res,isea3h_sp_pts@data[c("id","mu","var")]) %>%
        filter(centroid ==0) %>%
        clip_polygons_lonlat("id")

    g1 <- (LinePlotTheme() + geom_polygon(data=isea3h_clipped,aes(lon,lat,fill=mu,group=id),colour="light grey") + scale_fill_distiller(palette="Spectral",trans="reverse") + coord_map("mollweide")) %>%
        draw_world(inc_border=FALSE)

    mumin <- min(isea3h_clipped$mu)
    mumax <- max(isea3h_clipped$mu)

    g2 <- (LinePlotTheme() + geom_point(data=data.frame(AIRS_05_2003),aes(lon,lat,colour=pmin(pmax(co2avgret,mumin),mumax)),size=2) +
               scale_colour_distiller(palette="Spectral",trans="reverse",guide_legend(title="co2")) + coord_map("mollweide")) %>%
        draw_world(inc_border=TRUE)

    #g_all <- arrangeGrob(g1,g2,ncol=1)
})


