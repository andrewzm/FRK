context("other")

test_that("opts",{
    expect_is(opts_FRK$get,"function")
    expect_is(opts_FRK$set,"function")
    expect_is(opts_FRK$get("progress"),"logical")
    expect_is(opts_FRK$get("parallel"),"integer")
    expect_is(opts_FRK$get("verbose"),"logical")
})

test_that("coordinates",{
    library(sp)
    data <- data.frame(x = seq(0,10,by=0.01), y = 0, z= runif(1001))
    data_sp <- data
    coordinates(data_sp) <- ~x+y
    Grid1D_df <- auto_BAUs(manifold = real_line(),
                           cellsize = 1,
                           data=data_sp)
    test_coords <- coordinates(Grid1D_df)
    row.names(test_coords) <- NULL
    expect_equal(test_coords,as.matrix(data.frame(x=-2:12,y=0)))
})

test_that("coordnames_SpaceTime",{
    library(spacetime)

    sim_process <- expand.grid(x = seq(0.005,0.995,by=0.01),
                               y = seq(0.005,0.995,by=0.01),
                               t = seq(1,20,by = 1))
    sim_process$z <- 1

    time1 <- as.POSIXct("2015-09-01",tz="") + 3600*24*(sim_process$t-1)
    space1 <- sim_process[,c("x","y")]
    coordinates(space1) <- ~x+y
    STobj1 <- STIDF(space1,time1,data=sim_process)
    expect_equal(coordnames(STobj1),c("x","y","t"))

    time2 <- unique(time1)
    space2 <- unique(sim_process[,c("x","y")])
    coordinates(space2) <- ~x+y
    STobj2 <- STIDF(space2,time2,data=sim_process)
    expect_equal(coordnames(STobj2),c("x","y","t"))
})

test_that("Estimate observation error from variogram",{
    data(meuse)
    coordinates(meuse) = ~x+y # change into an sp object
    suppressWarnings(meuse <- .est_obs_error(meuse,variogram.formula = log(zinc)~1))
    expect_true("std" %in% names(meuse))
    expect_is(meuse,"SpatialPointsDataFrame")

})


test_that("Can convert SPDF to DF", {
    sphere_grid <- auto_BAUs(manifold=sphere(),
                             data=NULL,
                             cellsize=c(20,10),
                             type="grid")
    expect_is(sphere_grid,"SpatialPolygonsDataFrame")
    df <- SpatialPolygonsDataFrame_to_df(sphere_grid,vars = c("lon","lat"))
    expect_is(df,"data.frame")
    expect_equal(names(df),c("lon","lat","id"))
    expect_equal(length(unique(df$id)),length(sphere_grid))

    sphere_grid$test <- 1:length(sphere_grid)
    df <- SpatialPolygonsDataFrame_to_df(sphere_grid,vars = c("lon","lat","test"))
    expect_is(df,"data.frame")
    expect_equal(names(df),c("lon","lat","id","test"))
    expect_equal(length(unique(df$id)),length(sphere_grid))

})

test_that("Can convert DF to SP", {
    library(sp)
    opts_FRK$set("parallel",0L)
    df <- data.frame(id = c(rep(1,4),rep(2,4)),
                     x = c(0,1,0,0,2,3,2,2),
                     y=c(0,0,1,0,0,1,1,0))
    pols <- df_to_SpatialPolygons(df,"id",c("x","y"),CRS())
    expect_is(pols,"SpatialPolygons")
    expect_equal(length(pols),2)
    expect_equal(coordnames(pols),c("x","y"))

})

test_that("Options work", {
    opts <- new_opts_FRK()
    expect_is(opts,"list")
    expect_equal(names(opts),c("set","get"))
    expect_is(opts$set,"function")
    expect_is(opts$get,"function")
    opts$set("progress",1)
    expect_equal(opts$get("progress"),1)
})

test_that("Plotting works", {
    library(ggplot2)
    expect_true({draw_world(); TRUE})
})

test_that("Date sequencing works", {
    library("Hmisc")
    tspacing <- paste(1,"year")
    tstart <- as.POSIXct("1990-06-01 10:00:00 AEST")
    tend <- as.POSIXct("1993-06-01 10:00:00 AEST")
    tgrid <- seq(truncPOSIXt(tstart,"year"),
                 truncPOSIXt(tend,"year"),
                 by=tspacing)
    expect_is(tgrid,"POSIXct")
    expect_equal(length(tgrid),4L)
})

test_that("distR works", {
  x <- matrix(rnorm(20),10,2)
  y <- matrix(rnorm(20),10,2)
  D1 <- distR(x,x)
  D2 <- as.matrix(dist(x))
  expect_lt(sum(D1-D2), 1e-12)
})

test_that("formula no covariates works", {
    f1 <- y ~ x + 1
    f2 <- sin(y) ~ log(x) - 1
    f3 <- sin(cos(y)) ~ d + r

    f1B <- .formula_no_covars(f1)
    f2B <- .formula_no_covars(f2)
    f3B <- .formula_no_covars(f3)

    expect_equal(f1B,formula(y~1))
    expect_equal(f2B,formula(sin(y)~1))
    expect_equal(f3B,formula(sin(cos(y))~1))
})
