# Polygons from sp vignette as our observations
Sr1 = Polygon(cbind(c(2.03,4.06,4.06,2.03),c(2.03,2.03,4.06,4.06)))
Sr2 = Polygon(cbind(c(1.03,3.03,3.03,1.03),c(1.03,1.03,3.03,3.03)))
Srs1 = Polygons(list(Sr1), "s1")
Srs2 = Polygons(list(Sr2), "s2")
SpP = SpatialPolygons(list(Srs1,Srs2), 1:2)
attr = data.frame(z=c(5,2),std=c(1,1), row.names=c("s1", "s2"))
SrDf = SpatialPolygonsDataFrame(SpP, attr)
coordnames(SrDf) <- c("x","y")

# Now create the BAUs
BAUs <- auto_BAUs(manifold = plane(),
                  cellsize = c(0.3,0.3),
                  type = "grid",
                  data = SrDf,
                  convex=-0.3,
                  nonconvex_hull=FALSE)
idx1 <- which(cut(BAUs$x,c(2.03,4.06),labels = FALSE) &
                  cut(BAUs$y,c(2.03,4.06),labels = FALSE))
idx2 <- which(cut(BAUs$x,c(1.03,3.03),labels = FALSE) &
                  cut(BAUs$y,c(1.03,3.03),labels = FALSE))

test_that("observations with large support cover correct BAUs", {


    C <- BuildC(SrDf,BAUs)
    expect_equal(C$i_idx,c(idx1*0+1,idx2*0+2))
    expect_equal(C$j_idx,c(idx1,idx2))

    ## Visualise (not run)
    # BAUs@data$nn <- 0
    # BAUs@data$nn[C$j_idx[which(C$i_idx == 2)] ] <- 1

})

test_that("observations with large support average covariates correctly", {

    SrDf_updated <- map_data_to_BAUs(SrDf,BAUs,av_var=TRUE)

    ## Check that covariates are properly averaged and that overlapping obs. work
    obs1x <- mean(BAUs[["x"]][idx1])
    obs1y <- mean(BAUs[["y"]][idx1])


    obs2x <- mean(BAUs[["x"]][idx2])
    obs2y <- mean(BAUs[["y"]][idx2])

    expect_equal(obs1x,SrDf_updated@data$x[1])
    expect_equal(obs1y,SrDf_updated@data$y[1])
    expect_equal(obs2x,SrDf_updated@data$x[2])
    expect_equal(obs2y,SrDf_updated@data$y[2])

})
