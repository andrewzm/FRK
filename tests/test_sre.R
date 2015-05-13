library(sp)
library(ggplot2)
library(INLA)
library(dplyr)

print("Entering test_SRE.R")

## Planar SRE
test_that("can construct planar SRE", {
    data(meuse)
    meuse$std <- sqrt(0.05066)/5
    meuse$fs <- 1
    coordinates(meuse) = ~x+y # change into an sp object
    f <- log(zinc) ~ 1
    G <- auto_basis(m = plane(),data=meuse,nres = 2,prune=10,type = "bisquare")
    S <- SRE(f,meuse,G)
    g <- ggplot() + geom_point(data=data.frame(meuse),aes(x,y)); show_basis(g,S@basis) + coord_fixed()
})




## 1D SRE
test_that("can construct SRE on real line", {
    data(meuse)
    meuse$std <- sqrt(0.05066)
    meuse$fs <- 1
    meuse$y <- 0
    coordinates(meuse) = ~x+y # change into an sp object
    f <- log(zinc) ~ 1
    G <- auto_basis(m = real_line(),data=meuse,nres = 2,prune=4,type = "bisquare")
    S <- SRE(f,meuse,G)
    g <- ggplot() + geom_point(data=data.frame(meuse),aes(x,y)); show_basis(g,S@basis)
})

## Spherical SRE
test_that("can construct SRE on sphere", {
    N <- 1000
    random_data <- data.frame(lon = runif(n = N,min = -180,max=180),
                              lat = runif(n = N,min = -90, max = 90),
                              z = rnorm(N),
                              std = rnorm(N)/10,
                              fs = 1)
    coordinates(random_data) = ~lon + lat
    f <- z ~ lon + lat + 1
    G <- auto_basis(m = sphere(),data=random_data,nres = 3,prune=1,type = "bisquare")
    S <- SRE(f,random_data,G)
    g <- draw_world(EmptyTheme(),inc_border = T); show_basis(g,G) + coord_map("sinusoidal") + xlab("") + ylab("")
})

