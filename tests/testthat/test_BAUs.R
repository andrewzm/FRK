context("BAUs")

test_that("real_line_BAUs",{
    library(sp)
    data <- data.frame(x = seq(0,10,by=0.01), y = 0, z= runif(1001),std=0.5)
    coordinates(data) <- ~x+y
    Grid1D_df <- auto_BAUs(manifold = real_line(),
                           cellsize = 1,
                           data=data)
    expect_is(Grid1D_df,"SpatialPolygonsDataFrame")
    expect_equal(names(Grid1D_df),c("x","y"))
    expect_equal(mean(diff(Grid1D_df$x)),1)

    f <- z ~ 1
    binned_data <- map_data_to_BAUs(data,Grid1D_df,av_var=all.vars(f)[1])
    expect_is(binned_data,"SpatialPointsDataFrame")
    expect_true(nrow(binned_data) <= nrow(Grid1D_df))

    C <- BuildC(binned_data,Grid1D_df)
    expect_is(C,"list")
    expect_equal(names(C),c("i_idx","j_idx"))
    expect_equal(length(C$i_idx),nrow(binned_data))
    expect_equal(length(C$j_idx),nrow(binned_data))


})


test_that("plane_BAUs",{
    library(sp)
    data <- data.frame(x = rnorm(5),y=rnorm(5),z = rnorm(5),std=1)
    coordinates(data) <- ~x+y
    Grid2D <- auto_BAUs(manifold = plane(),
                           type="grid",
                           cellsize = 0.5,
                           data=data)
    expect_is(Grid2D,"SpatialPolygonsDataFrame")
    expect_equal(names(Grid2D),c("x","y"))


    f <- z ~ 1
    binned_data <- map_data_to_BAUs(data,Grid2D,av_var=all.vars(f)[1])
    expect_is(binned_data,"SpatialPointsDataFrame")
    expect_true(nrow(binned_data) <= nrow(Grid2D))

    C <- BuildC(binned_data,Grid2D)
    expect_is(C,"list")
    expect_equal(names(C),c("i_idx","j_idx"))
    expect_equal(length(C$i_idx),nrow(binned_data))
    expect_equal(length(C$j_idx),nrow(binned_data))
})


test_that("sphere_BAUs",{
    isea3h_1 <- auto_BAUs(manifold=sphere(),
                                 res=1,
                                 data=NULL,
                                 type="hex")
    expect_is(isea3h_1,"SpatialPolygonsDataFrame")
    expect_equal(nrow(isea3h_1@data),23)
    expect_equal(names(isea3h_1@data),c("id","lon","lat"))
    expect_true(grepl("+proj=longlat",proj4string(isea3h_1)))

    sphere_grid <- auto_BAUs(manifold=sphere(),
                             data=NULL,
                             cellsize=c(20,10),
                             type="grid")
    expect_is(sphere_grid,"SpatialPolygonsDataFrame")
    expect_equal(nrow(sphere_grid@data),324)
    expect_equal(names(sphere_grid@data),c("lon","lat","id"))
    expect_true(grepl("+proj=longlat",proj4string(sphere_grid)))

})


test_that("SpaceTime_BAUs",{
    library(sp)
    library(spacetime)
    sim_process <- expand.grid(x = seq(0.005,0.995,by=0.1),
                               y = seq(0.005,0.995,by=0.1),
                               t = seq(1,5,by = 1),
                               std = 0.5)
    sim_process$z <- 1

    time1 <- as.POSIXct("2015-09-01",tz="") + 3600*24*(sim_process$t-1)
    space1 <- sim_process[,c("x","y")]
    coordinates(space1) <- ~x+y
    STobj1 <- STIDF(space1,time1,data=sim_process)

    time_grid <- auto_BAUs(timeline(),cellsize = 1,d = STobj1,tunit="days")
    expect_is(time_grid,"POSIXt")

    space_time_grid <- auto_BAUs(STplane(),cellsize = c(0.4,0.4,1),type="hex",data = STobj1,tunit="days")
    expect_is(space_time_grid,"STFDF")

    f <- z ~ 1
    binned_data <- map_data_to_BAUs(STobj1,space_time_grid,av_var=all.vars(f)[1])
    expect_is(binned_data,"STIDF")

    C <- BuildC(binned_data,space_time_grid)
    expect_is(C,"list")
    expect_equal(names(C),c("i_idx","j_idx"))
    expect_equal(length(C$i_idx),as.numeric(nrow(binned_data)))
    expect_equal(length(C$j_idx),as.numeric(nrow(binned_data)))


})

