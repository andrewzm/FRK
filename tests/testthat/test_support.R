test_that("observations with large support are average covariates correctly", {

    # Polygons from sp vignette as our observations
    Sr1 = Polygon(cbind(c(2.03,4.06,4.06,2.03),c(2.03,2.03,4.06,4.06)))
    Sr2 = Polygon(cbind(c(1.03,3.03,3.03,1.03),c(1.03,1.03,3.03,3.03)))
    Srs1 = Polygons(list(Sr1), "s1")
    Srs2 = Polygons(list(Sr2), "s2")
    SpP = SpatialPolygons(list(Srs1,Srs2), 1:2)
    attr = data.frame(z=c(5,2),std=c(1,1), row.names=c("s1", "s2"))
    SrDf = SpatialPolygonsDataFrame(SpP, attr)

    # Now create the BAUs
    BAUs <- auto_BAUs(manifold = plane(),
                      cellsize = c(0.3,0.3),
                      type = "grid",
                      data = SrDf,
                      convex=-0.3)

    SrDf_updated <- map_data_to_BAUs(SrDf,BAUs,av_var=TRUE)

    BAU_as_points <- SpatialPointsDataFrame(coordinates(BAUs),BAUs@data)

    ## Check that covariates are properly averaged and that overlapping obs. work
    idx <- which(cut(BAUs$x,c(2.03,4.06),labels = FALSE) &
                     cut(BAUs$y,c(2.03,4.06),labels = FALSE))
    obs1x <- mean(BAUs[["x"]][idx])
    obs1y <- mean(BAUs[["y"]][idx])

    idx <- which(cut(BAUs$x,c(1.03,3.03),labels = FALSE) &
                     cut(BAUs$y,c(1.03,3.03),labels = FALSE))
    obs2x <- mean(BAUs[["x"]][idx])
    obs2y <- mean(BAUs[["y"]][idx])

    expect_equal(obs1x,SrDf_updated@data$x[1])
    expect_equal(obs1y,SrDf_updated@data$y[1])
    expect_equal(obs2x,SrDf_updated@data$x[2])
    expect_equal(obs2y,SrDf_updated@data$y[2])
})
