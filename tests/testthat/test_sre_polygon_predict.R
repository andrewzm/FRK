test_that("can predict over polygons in plane", {

    ## Get data
    library("sp")
    library("dplyr")

    data("meuse")
    meuse$std <- sqrt(0.05066)
    coordinates(meuse) = ~x+y # change into an sp object
    data("meuse.grid")
    gridded(meuse.grid) = ~x + y
    HexPts <- spsample(meuse.grid, type = "hexagonal", cellsize = 200)
    HexPols <- HexPoints2SpatialPolygons(HexPts)
    HexPts_df <- SpatialPointsDataFrame(HexPts,data.frame(fs=rep(1,length(HexPts))))
    HexPols_df <- SpatialPolygonsDataFrame(HexPols,
                                           over(HexPols,meuse.grid))

    G <- auto_basis(manifold = plane(),data=meuse,nres = 2,prune=10,type = "Gaussian")
    expect_is(G, "Basis")

    S1 <- eval_basis(G,HexPts_df)
    S2 <- eval_basis(G,HexPols_df)
    expect_is(S1, "Matrix")
    expect_is(S2, "Matrix")
    expect_identical(dim(S1), c(length(HexPols), nbasis(G)))
    expect_identical(dim(S2), c(length(HexPols), nbasis(G)))
    #plot(as.numeric(S1))
    #lines(as.numeric(S2),col='red')

    HexPols_df <- auto_BAUs(manifold = plane(),
                            cellsize = c(400,400),
                            type = "grid",
                            data = meuse,
                            convex=-0.05,
                            nonconvex_hull=FALSE)
    HexPols_df$fs <- 1
    expect_is(HexPols_df, "SpatialPixelsDataFrame")

    f <- log(zinc) ~ 1
    S <- SRE(f,data = list(meuse),
             basis = G,
             BAUs = HexPols_df,
             est_error = FALSE)
    S <- SRE.fit(S, n_EM = 10, print_lik=F)
    expect_is(S, "SRE")

    HexPols_df <- predict(S)
    expect_is(HexPols_df, "SpatialPixelsDataFrame")
    expect_true(all(c("var","mu","sd") %in% names(HexPols_df)))

    ## Try with sigma2fs = 0
    S@sigma2fshat <- 0
    HexPols_df2 <- predict(S)
    expect_is(HexPols_df2, "SpatialPixelsDataFrame")
    expect_true(all(c("var","mu","sd") %in% names(HexPols_df2)))
    #spplot(HexPols_df,"mu")
})

