test_that("SRE 1D works",{

    library(dplyr)
    library(sp)
    ### Generate process and data
    sim_process <- data.frame(x = seq(0.005,0.995,by=0.01)) %>%
        mutate(y=0,proc = sin(x*10) + 0.3*rnorm(length(x)))
    sim_data <- sample_n(sim_process,50) %>%
        mutate(z = proc + 0.1*rnorm(length(x)), std = 0.1)
    coordinates(sim_data) = ~x + y# change into an sp object
    grid_BAUs <- auto_BAUs(manifold=real_line(),data=sim_data,cellsize = c(0.01),type="grid")
    grid_BAUs$fs = 1

    ### Set up SRE model
    G <- auto_basis(m = real_line(),
                    data=sim_data,
                    nres = 2,
                    regular = 6,
                    type = "bisquare",
                    subsamp = 20000)
    expect_is(G,"Basis")

    f <- z ~ 1
    S <- SRE(f,list(sim_data),G,
             grid_BAUs,
             est_error = FALSE)

    expect_is(S,"SRE")

    ### Fit with 5 EM iterations so as not to take too much time
    S <- SRE.fit(S,n_EM = 3,tol = 1e-5,print_lik=TRUE)
    expect_is(S,"SRE")

    ### Predict over BAUs
    grid_BAUs <- SRE.predict(S,pred_locs = grid_BAUs,use_centroid = TRUE)
    expect_is(grid_BAUs,"SpatialPolygonsDataFrame")

    ### summary works?
    expect_true({summary(S);TRUE})
})

test_that("SRE 2D plane works",{

    library(dplyr)
    library(sp)
    ### Generate process and data
    sim_process <- expand.grid(x = seq(0,1,by=0.2),y = seq(0,1,by=0.2)) %>%
        mutate(proc = sin(x*10) + 0.3*rnorm(length(x)))
    sim_data <- sample_n(sim_process,10) %>%
        mutate(z = proc + 0.1*rnorm(length(x)), std = 0.1)
    coordinates(sim_data) = ~x + y# change into an sp object
    grid_BAUs <- auto_BAUs(manifold=plane(),data=sim_data,cellsize = c(0.4),type="grid")
    grid_BAUs$fs = 1

    ### Set up SRE model
    G <- auto_basis(m = plane(),
                    data=sim_data,
                    nres = 1,
                    regular = 2,
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
    grid_BAUs <- SRE.predict(S,pred_locs = grid_BAUs,use_centroid = TRUE)
    expect_is(grid_BAUs,"SpatialPolygonsDataFrame")

    ### summary works?
    expect_true({summary(S);TRUE})
})


test_that("SRE sphere works",{

    library(dplyr)
    library(sp)
    ### Generate process and data
    sim_process <- expand.grid(lon = seq(-180,180,by=10),lat = seq(-90,90,by=10)) %>%
        mutate(proc = sin(lon*10) + 0.3*rnorm(length(lon)))
    sim_data <- sample_n(sim_process,100) %>%
        mutate(z = proc + 0.1*rnorm(length(lon)), std = 0.1)
    coordinates(sim_data) = ~lon + lat# change into an sp object
    proj4string(sim_data)=CRS("+proj=longlat +ellps=sphere")
    grid_BAUs <- auto_BAUs(manifold=sphere(),data=sim_data,res = 2,type="hex")
    grid_BAUs$fs = 1

    ### Set up SRE model
    G <- auto_basis(m = sphere(),
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
    S <- SRE.fit(S,n_EM = 3,tol = 1e-5,print_lik=TRUE)
    expect_is(S,"SRE")

    ### Predict over BAUs
    grid_BAUs <- SRE.predict(S,pred_locs = grid_BAUs,use_centroid = TRUE)
    grid_BAUs <- SRE.predict(S,pred_locs = grid_BAUs,use_centroid = FALSE)
    expect_is(grid_BAUs,"SpatialPolygonsDataFrame")

    ### summary works?
    expect_true({summary(S);TRUE})
})

test_that("SRE space-time sphere works",{

    library(dplyr)
    library(sp)
    library(spacetime)
    ### Generate process and data
    sim_process <- expand.grid(lon = seq(-180,180,by=20),lat = seq(-90,90,by=20),t=1:5) %>%
        mutate(proc = sin(lon*10) + 0.3*rnorm(length(lon)))
    sim_data <- sample_n(sim_process,100) %>%
        mutate(z = proc + 0.1*rnorm(length(lon)), std = 0.1)

    time <- as.POSIXct("2003-05-01",tz="") + 3600*24*(sim_data$t-1)
    space <- sim_data[,c("lon","lat")]
    coordinates(space) = ~lon+lat # change into an sp object
    proj4string(space)=CRS("+proj=longlat +ellps=sphere")
    STobj <- STIDF(space,time,data=sim_data)

    grid_BAUs <- auto_BAUs(manifold=STsphere(),res=1,data=STobj,type="hex")
    grid_BAUs$fs = 1


    ### Set up SRE model
    G_spatial <- auto_basis(m = sphere(),
                            data=as(STobj,"Spatial"),
                            nres = 1,
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
    S <- SRE.fit(S,n_EM = 3,tol = 1e-5,print_lik=TRUE)
    expect_is(S,"SRE")

    ### Predict over BAUs
    grid_BAUs <- SRE.predict(S,pred_locs = grid_BAUs,use_centroid = TRUE)
    expect_is(grid_BAUs,"STFDF")

    ### summary works?
    expect_true({summary(S);TRUE})
})


