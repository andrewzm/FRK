test_that("SRE 1D works",{

    library(dplyr)
    library(sp)
    library(Matrix)
    ### Generate process and data
    sim_process <- data.frame(x = seq(0.005,0.995,by=0.01)) %>%
        mutate(y=0,proc = sin(x*10) + 0.3*rnorm(length(x)))
    sim_data <- sample_n(sim_process,50) %>%
        mutate(z = proc + 0.1*rnorm(length(x)), std = 0.1)
    coordinates(sim_data) = ~x + y# change into an sp object
    grid_BAUs <- auto_BAUs(manifold=real_line(),
                           type="grid",
                           data=sim_data,
                           cellsize = c(0.01))
    grid_BAUs$fs = 1

    ### Set up SRE model
    G <- auto_basis(manifold = real_line(),
                    data=sim_data,
                    nres = 2,
                    regular = 6,
                    type = "bisquare",
                    subsamp = 20000)
    expect_is(G,"Basis")

    c_res <- count_res(G)
    expect_is(c_res,"data.frame")
    expect_equal(c_res$res,c(1,2))

    f <- z ~ 1
    S <- SRE(f,list(sim_data),G,
             grid_BAUs,
             est_error = FALSE)

    c_res2 <- count_res(S)
    expect_identical(c_res,c_res2)
    expect_is(S,"SRE")

    ### Fit with 5 EM iterations so as not to take too much time
    S <- SRE.fit(S,n_EM = 3,tol = 1e-5,print_lik=FALSE)
    expect_is(S,"SRE")
    expect_is(info_fit(S),"list")

    ### Predict over BAUs using both modes
    grid_BAUs <- predict(S, covariances = TRUE)
    expect_is(grid_BAUs,"list")
    expect_equal(as.numeric(grid_BAUs$newdata$var), unname(diag(grid_BAUs$Cov)))
    grid_BAUs <- predict(S, obs_fs = FALSE, covariances = TRUE)
    expect_is(grid_BAUs,"list")
    expect_equal(as.numeric(grid_BAUs$newdata$var), unname(diag(grid_BAUs$Cov)))
    expect_equal(dim(grid_BAUs$Cov),rep(length(S@BAUs),2))

    ### Check covariances option
    grid_BAUs <- predict(S)
    expect_is(grid_BAUs,"SpatialPixelsDataFrame")
    grid_BAUs <- predict(S,obs_fs = FALSE)
    expect_is(grid_BAUs,"SpatialPixelsDataFrame")
    print(coef(S))

    ### Check alphahat
    alphahat <- coef(S)
    expect_is(alphahat, "numeric")
    expect_equal(length(alphahat), 1L)
    expect_equal(names(alphahat), "Intercept")

    ### summary/print/show works
    expect_true({summary(S);TRUE})
    expect_true({print(S);TRUE})
    expect_true({show(S);TRUE})
})

test_that("SRE 2D plane works",{

    library(dplyr)
    library(sp)
    set.seed(1)
    ### Generate process and data
    sim_process <- expand.grid(x = seq(0,1,by=0.2),y = seq(0,1,by=0.2)) %>%
        mutate(proc = sin(x*10) + 0.3*rnorm(length(x)))
    sim_data <- sample_n(sim_process,10) %>%
        mutate(z = proc + 0.1*rnorm(length(x)), std = 0.1)
    coordinates(sim_data) = ~x + y# change into an sp object
    grid_BAUs <- auto_BAUs(manifold=plane(),
                           type="grid",
                           data=sim_data,
                           cellsize = c(0.4),
                           nonconvex_hull=FALSE)
    grid_BAUs$fs = 1

    ### Set up SRE model
    G <- auto_basis(manifold = plane(),
                    data=sim_data,
                    nres = 1,
                    regular = 1,
                    type = "bisquare",
                    subsamp = 20000)
    expect_is(G,"Basis")

    f <- z ~ 1
    S <- SRE(f,list(sim_data),G,
             grid_BAUs,
             est_error = FALSE)

    expect_is(S,"SRE")

    ### Fit with 5 EM iterations so as not to take too much time
    S <- SRE.fit(S,n_EM = 3,tol = 1e-5,print_lik=FALSE)
    expect_is(S,"SRE")

    ### Predict over BAUs
    grid_BAUs <- predict(S)
    expect_is(grid_BAUs,"SpatialPixelsDataFrame")
    grid_BAUs <- predict(S,obs_fs = FALSE)
    expect_is(grid_BAUs,"SpatialPixelsDataFrame")

    ### summary works?
    expect_true({summary(S);TRUE})
})


test_that("SRE sphere works",{

    library(dplyr)
    library(sp)
    ### Generate process and data
    sim_process <- expand.grid(lon = seq(-100,100,by=10),lat = seq(-50,50,by=10)) %>%
        mutate(proc = sin(lon*10) + 0.3*rnorm(length(lon)))
    sim_data <- sample_n(sim_process,100) %>%
        mutate(z = proc + 0.1*rnorm(length(lon)), std = 0.1)
    coordinates(sim_data) = ~lon + lat# change into an sp object
    slot(sim_data, "proj4string") = CRS("+proj=longlat +ellps=sphere")
    grid_BAUs <- auto_BAUs(manifold = sphere(),
                           type = "hex",
                           data = sim_data,
                           isea3h_res = 2)
    grid_BAUs$fs = 1

    ### Set up SRE model
    G <- auto_basis(manifold = sphere(),
                    data=sim_data,
                    nres = 1,
                    type = "bisquare",
                    subsamp = 20000)
    expect_is(G,"Basis")

    f <- z ~ 1
    S <- SRE(f,list(sim_data),G,
             grid_BAUs,
             est_error = FALSE)

    expect_is(S,"SRE")

    ### Fit with 5 EM iterations so as not to take too much time
    S <- SRE.fit(S,n_EM = 3,tol = 1e-5,print_lik=FALSE)
    expect_is(S,"SRE")

    ### Predict over BAUs
    grid_BAUs <- predict(S)
    expect_is(grid_BAUs,"SpatialPolygonsDataFrame")
    grid_BAUs <- predict(S,obs_fs = FALSE)
    expect_is(grid_BAUs,"SpatialPolygonsDataFrame")

    ### summary works?
    expect_true({summary(S);TRUE})
})

test_that("SRE space-time sphere works",{

    library(dplyr)
    library(sp)
    library(spacetime)
    ### Generate process and data
    set.seed(1)
    sim_process <- expand.grid(lon = seq(-100,100,by=20),lat = seq(-40,40,by=20),t=1:5) %>%
        mutate(proc = sin(lon*10) + 0.3*rnorm(length(lon)))
    sim_data <- sample_n(sim_process,100) %>%
        mutate(z = proc + 0.1*rnorm(length(lon)), std = 0.1)

    time <- as.POSIXct("2003-05-01",tz="") + 3600*24*(sim_data$t-1)
    space <- sim_data[,c("lon","lat")]
    coordinates(space) = ~lon+lat # change into an sp object
    proj4string(space)=CRS("+proj=longlat +ellps=sphere")
    sim_data$t <- NULL
    STobj <- STIDF(space,time,data=sim_data)

    grid_BAUs <- auto_BAUs(manifold=STsphere(),
                           type="hex",
                           isea3h_res=1,
                           data=STobj)
    grid_BAUs$fs = 1


    ### Set up SRE model
    G_spatial <- auto_basis(manifold = sphere(),
                            data=as(STobj,"Spatial"),
                            nres = 2,
                            type = "bisquare",
                            subsamp = 20000)

    G_temporal <- local_basis(manifold=real_line(),loc = matrix(c(1,3)),scale = rep(1,2))
    G <- TensorP(G_spatial,G_temporal)
    expect_is(G,"TensorP_Basis")

    f <- z ~ 1
    S <- SRE(f,list(STobj),G,
             grid_BAUs,
             est_error = FALSE)

    expect_is(S,"SRE")

    ### Fit with 5 EM iterations so as not to take too much time
    S <- SRE.fit(S,n_EM = 3,tol = 1e-5,print_lik=FALSE)
    expect_is(S,"SRE")

    ### Predict over BAUs
    grid_BAUs <- predict(S)
    expect_is(grid_BAUs,"STFDF")
    grid_BAUs <- predict(S,obs_fs = FALSE)
    expect_is(grid_BAUs,"STFDF")

    ### summary works?
    expect_true({summary(S);TRUE})
    
    ## Test with multiple ST objects with differing amounts of data
    point_to_remove <- 2
    updated_data <- STobj@data[-point_to_remove, ]
    updated_coords <- STobj@sp[-point_to_remove, ]
    updated_time <- STobj@time[-point_to_remove]
    STobj2 <- STIDF(updated_coords, updated_time, data = updated_data)
    S <- SRE(f,list(STobj, STobj2),G,
             grid_BAUs,
             est_error = FALSE)
    expect_is(S,"SRE")
    S <- SRE.fit(S,n_EM = 3,tol = 1e-5,print_lik=FALSE)
    expect_is(S,"SRE")
    grid_BAUs <- predict(S)
    expect_is(grid_BAUs,"STFDF")
    grid_BAUs <- predict(S,obs_fs = FALSE)
    expect_is(grid_BAUs,"STFDF")
    expect_true({summary(S);TRUE})
})


