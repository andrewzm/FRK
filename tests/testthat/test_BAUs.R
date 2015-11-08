context("BAUs")

test_that("real_line_BAUs",{
    library(sp)
    data <- data.frame(x = seq(0,10,by=0.01), y = 0, z= runif(1001))
    coordinates(data) <- ~x+y
    Grid1D_df <- auto_BAUs(manifold = real_line(),
                           cellsize = 1,
                           data=data)
    expect_is(Grid1D_df,"SpatialPolygonsDataFrame")
    expect_equal(names(Grid1D_df),c("x","y"))
    expect_equal(mean(diff(Grid1D_df$x)),1)
})

test_that("sphere_BAUs",{
    isea3h_1 <- auto_BAUs(manifold=sphere(),
                                 res=1,
                                 data=NULL,
                                 type="hex")
    expect_is(isea3h_1,"SpatialPolygonsDataFrame")
    expect_equal(nrow(isea3h_1@data),23)
    expect_equal(names(isea3h_1@data),c("id","lon","lat"))
    expect_equal(proj4string(isea3h_1),"+proj=longlat")

    sphere_grid <- auto_BAUs(manifold=sphere(),
                             data=NULL,
                             cellsize=c(20,10),
                             type="grid")
    expect_is(sphere_grid,"SpatialPolygonsDataFrame")
    expect_equal(nrow(sphere_grid@data),324)
    expect_equal(names(sphere_grid@data),c("lon","lat","id"))
    expect_equal(proj4string(isea3h_1),"+proj=longlat")

})
