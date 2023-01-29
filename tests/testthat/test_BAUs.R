context("BAUs")

test_that("real_line_BAUs",{
    library(sp)
    data <- data.frame(x = seq(0,10,by=0.01), y = 0, z= runif(1001),std=0.5)
    coordinates(data) <- ~x+y
    Grid1D_df <- auto_BAUs(manifold = real_line(),
                           cellsize = 1,
                           data=data)
    expect_is(Grid1D_df,"SpatialPixelsDataFrame")
    expect_equal(names(Grid1D_df),c("x","y"))
    expect_equal(mean(diff(Grid1D_df$x)),1)

    f <- z ~ 1
    binned_data1 <- map_data_to_BAUs(data,Grid1D_df,average_in_BAU = TRUE)
    binned_data2 <- map_data_to_BAUs(data,Grid1D_df,average_in_BAU = FALSE)
    expect_is(binned_data1,"SpatialPointsDataFrame")
    expect_is(binned_data2,"SpatialPointsDataFrame")
    expect_true(nrow(binned_data1) <= nrow(Grid1D_df))

    C1 <- BuildC(binned_data1,Grid1D_df)
    C2 <- BuildC(binned_data2,Grid1D_df)
    expect_is(C1,"list")
    expect_equal(names(C1),c("i_idx","j_idx", "x_idx"))
    expect_equal(length(C1$i_idx),nrow(binned_data1))
    expect_equal(length(C1$j_idx),nrow(binned_data1))

})

test_that("plane_BAUs",{
    library(sp)
    set.seed(1)
    data <- data.frame(x = rnorm(5),y=rnorm(5),z = rnorm(5),std=1)
    coordinates(data) <- ~x+y
    if(require("INLA") & require("rgdal", quietly = TRUE)) {
        Grid2D <- auto_BAUs(manifold = plane(),
                            type="grid",
                            cellsize = 0.5,
                            data=data,
                            nonconvex_hull = TRUE)
        expect_is(Grid2D,"SpatialPixelsDataFrame")
        expect_equal(names(Grid2D),c("x","y"))
    }


    ## Now without INLA
    Grid2D <- auto_BAUs(manifold = plane(),
                        type="grid",
                        cellsize = 0.5,
                        data=data,
                        nonconvex_hull = FALSE)
    expect_is(Grid2D,"SpatialPixelsDataFrame")
    expect_equal(names(Grid2D),c("x","y"))


    f <- z ~ 1
    binned_data <- map_data_to_BAUs(data,Grid2D)
    expect_is(binned_data,"SpatialPointsDataFrame")
    expect_true(nrow(binned_data) <= nrow(Grid2D))

    C <- BuildC(binned_data,Grid2D)
    expect_is(C,"list")
    expect_equal(names(C),c("i_idx","j_idx", "x_idx"))
    expect_equal(length(C$i_idx),nrow(binned_data))
    expect_equal(length(C$j_idx),nrow(binned_data))

    ## Limited 2D grid
    Grid2D_limited <- auto_BAUs(manifold = plane(),
                        type="grid",
                        cellsize = 0.5,
                        data=data,
                        nonconvex_hull = FALSE,
                        xlims=c(-2,2),
                        ylims=c(-2,2))
    expect_is(Grid2D_limited,"SpatialPixelsDataFrame")
    expect_equal(names(Grid2D_limited),c("x","y"))
    expect_equal(min(Grid2D_limited@data[,1]),-2)
    expect_equal(max(Grid2D_limited@data[,1]),2)
    expect_equal(min(Grid2D_limited@data[,2]),-2)
    expect_equal(max(Grid2D_limited@data[,2]),2)
})


test_that("sphere_BAUs",{
    isea3h_1 <- auto_BAUs(manifold=sphere(),
                          type="hex",
                          isea3h_res=1,
                          data=NULL)
    expect_is(isea3h_1,"SpatialPolygonsDataFrame")
    expect_equal(nrow(isea3h_1@data),39)
    expect_equal(names(isea3h_1@data),c("id","lon","lat"))
    expect_true(grepl("+proj=longlat",.rawproj4string(isea3h_1)))

    sphere_grid <- auto_BAUs(manifold=sphere(),
                             type="grid",
                             data=NULL,
                             cellsize=c(20,10))
    expect_is(sphere_grid,"SpatialPolygonsDataFrame")
    expect_equal(nrow(sphere_grid@data),324)
    expect_equal(names(sphere_grid@data),c("lon","lat"))
    expect_true(grepl("+proj=longlat",.rawproj4string(sphere_grid)))

    sphere_grid_limited <- auto_BAUs(manifold=sphere(),
                             type="grid",
                             data=NULL,
                             cellsize=c(20,10),
                             xlims=c(-100,120),
                             ylims=c(-80,70))
    expect_is(sphere_grid_limited,"SpatialPolygonsDataFrame")
    expect_equal(nrow(sphere_grid_limited@data),165)
    expect_equal(names(sphere_grid_limited@data),c("lon","lat"))
    expect_true(grepl("+proj=longlat",.rawproj4string(sphere_grid_limited)))
    expect_equal(min(sphere_grid_limited@data[,1]),-90)
    expect_equal(max(sphere_grid_limited@data[,1]),110)
    expect_equal(min(sphere_grid_limited@data[,2]),-75)
    expect_equal(max(sphere_grid_limited@data[,2]),65)

})

test_that("sphere_BAUs_subset_BAUs",{

    set.seed(1)
    df <- data.frame(lon = runif(n = 1000, min = 120, max = 160),
                     lat = runif(n = 1000, min = 57, max = 88))
    coordinates(df) <- c("lon", "lat")
    slot(df, "proj4string") <- CRS('+proj=longlat +ellps=sphere')

    isea3h_1 <- auto_BAUs(manifold = sphere(),
                          type = "hex",
                          isea3h_res = 5,
                          data = df)

    sf::sf_use_s2(FALSE)
    data_in_BAUs <- sf::st_contains(as(isea3h_1, "sf"), as(df, "sf"))
    expect_equal(all(colSums(as.matrix(data_in_BAUs)) == 1), TRUE)
    #plot(isea3h_1, col = "red")
    #plot(df, add = TRUE)
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

    time_grid <- auto_BAUs(real_line(),
                           cellsize = 1,
                           d = as.Date(time(STobj1)),
                           tunit="days")
    expect_is(time_grid,"POSIXct")

    if(require("INLA") & require("rgdal", quietly = TRUE)) {
        space_time_grid <- auto_BAUs(STplane(),
                                     type="hex",
                                     cellsize = c(0.1,0.1,1),
                                     data = STobj1,
                                     tunit="days",
                                     convex= -0.2,
                                     nonconvex_hull = TRUE)
        expect_is(space_time_grid,"STFDF")
        expect_is(time(space_time_grid),"POSIXct")
    }

    space_time_grid <- auto_BAUs(STplane(),
                                 type="hex",
                                 cellsize = c(0.1,0.1,1),
                                 data = STobj1,
                                 tunit="days",
                                 convex= -0.2,
                                 nonconvex_hull = FALSE)
    expect_is(space_time_grid,"STFDF")
    expect_is(time(space_time_grid),"POSIXct")

    STobj2 <- space_time_grid[1:5,1:3] # mock space-time STFDF data
    STobj2$z <- 1

    f <- z ~ 1
    binned_data1 <- FRK:::map_data_to_BAUs(STobj1,
                                           space_time_grid,
                                           average_in_BAU = TRUE)
    binned_data2 <- FRK:::map_data_to_BAUs(STobj1,
                                           space_time_grid,
                                           average_in_BAU = FALSE)
    expect_true(ncol(binned_data2) >= ncol(binned_data1))
    expect_is(binned_data1,"STIDF")
    expect_is(binned_data2,"STIDF")

    C1 <- BuildC(binned_data1,space_time_grid)
    C2 <- BuildC(binned_data2,space_time_grid)
    expect_is(C1,"list")
    expect_is(C2,"list")
    expect_equal(names(C1),c("i_idx","j_idx", "x_idx"))
    expect_equal(names(C2),c("i_idx","j_idx", "x_idx"))
    expect_equal(length(C1$i_idx),as.numeric(nrow(binned_data1)))
    expect_equal(length(C1$j_idx),as.numeric(nrow(binned_data1)))

    ## Now do the same but with a slightly shifted time
    STobj3 <- STIDF(space1,time1 + 4000,data=sim_process)
    STobj3$z <- 1

    f <- z ~ 1
    binned_data3 <- FRK:::map_data_to_BAUs(STobj3,
                                           space_time_grid,
                                           average_in_BAU = TRUE)
    binned_data4 <- FRK:::map_data_to_BAUs(STobj3,
                                           space_time_grid,
                                           average_in_BAU = FALSE)
    expect_true(ncol(binned_data4) >= ncol(binned_data3))
    expect_is(binned_data3,"STIDF")
    expect_is(binned_data4,"STIDF")

    space_time_grid2 <- auto_BAUs(STplane(),
                                 type="hex",
                                 cellsize = c(0.1,0.1,1),
                                 data = STobj3,
                                 tunit="days",
                                 convex= -0.2,
                                 nonconvex_hull = FALSE)

    expect_equal(attr(space_time_grid@time,"tzone"),attr(STobj1@time,"tzone"))
})


test_that("Point from BAUs works",{
    library(sp)
    dat <- data.frame(x = rnorm(100),
                      y = rnorm(100))
    coordinates(dat) <- ~x+y
    BAUs <- BAUs_from_points(dat)

    expect_is(BAUs,"SpatialPolygonsDataFrame")
    expect_equal(length(BAUs),100)

    dat$z <- rnorm(100)
    BAUs <- BAUs_from_points(dat)
    expect_is(BAUs,"SpatialPolygonsDataFrame")
    expect_equal(length(BAUs),100)
})
