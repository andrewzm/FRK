test_that("can predict over polygons in plane", {

    ## Get data
    library(sp)
    library(dplyr)

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
                                            over(HexPols,meuse) %>%
                                                select(zinc)) %>%
                    subset(!is.na(zinc))


    G <- auto_basis(m = plane(),data=meuse,nres = 2,prune=10,type = "Gaussian")

    S1 <- eval_basis(G,HexPts_df)
    S2 <- eval_basis(G,HexPols_df)
    plot(as.numeric(S1))
    lines(as.numeric(S2),col='red')

    HexPols_df <- auto_BAUs(manifold = plane(),
                            cellsize = c(400,400),
                            type = "grid",
                            data = meuse,
                            convex=-0.05)
    HexPols_df$fs <- 1
    f <- log(zinc) ~ 1
    S <- SRE(f,data = list(meuse),basis = G,BAUs = HexPols_df,est_error = FALSE)
    S <- SRE.fit(S,n_EM = 10,print_lik=T)
    HexPols_df <- SRE.predict(S,HexPols_df,use_centroid = FALSE)
    spplot(HexPols_df,"mu")
})

