context("other")
test_that("opts",{
    expect_is(opts_FRK$get,"function")
    expect_is(opts_FRK$set,"function")
    expect_is(opts_FRK$get("progress"),"logical")
    expect_is(opts_FRK$get("parallel"),"integer")
    expect_is(opts_FRK$get("verbose"),"logical")
    expect_is(opts_FRK$get("Rhipe"),"logical")
})

test_that("coordinates",{
    library(sp)
    data <- data.frame(x = seq(0,10,by=0.01), y = 0, z= runif(1001))
    data_sp <- data
    coordinates(data_sp) <- ~x+y
    Grid1D_df <- auto_BAUs(manifold = real_line(),
                           cellsize = 1,
                           data=data_sp)
    expect_identical(coordinates(Grid1D_df),as.matrix(data.frame(x=-2:12,y=0)))
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
    STobj2 <- STIDF(space2,time2,data=sim_data)
    expect_equal(coordnames(STobj2),c("x","y","t"))
})

test_that("Estimate observation error from variogram",{
    data(meuse)
    meuse$Nobs <- 1
    coordinates(meuse) = ~x+y # change into an sp object
    suppressWarnings(meuse <- est_obs_error(meuse,variogram.formula = log(zinc)~1))
    expect_true("std" %in% names(meuse))
    expect_is(meuse,"SpatialPointsDataFrame")

})
