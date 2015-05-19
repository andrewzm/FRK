library(sp)
library(ggplot2)
library(INLA)
library(dplyr)

print("Entering test_sre_fit.R")

## Planar SRE
test_that("can fit planar SRE", {
    data(meuse)
    meuse$std <- sqrt(0.05066)/5
    meuse$fs <- 1
    coordinates(meuse) = ~x+y # change into an sp object
    f <- log(zinc) ~ 1
    G <- auto_basis(m = plane(),data=meuse,nres = 2,prune=10,type = "bisquare")
    S <- SRE(f,meuse,G)
    #g <- ggplot() + geom_point(data=data.frame(meuse),aes(x,y)); show_basis(g,S@basis) + coord_fixed()
})

