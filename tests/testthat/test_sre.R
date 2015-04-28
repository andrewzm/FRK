library(sp)
library(ggplot2)
library(INLA)
library(dplyr)

data(meuse)
meuse$std <- sqrt(0.05066)
coordinates(meuse) = ~x+y # change into an sp object
f <- log(zinc) ~ 1

G <- auto_basis(m = plane(),data=meuse,nres = 2,type = "bisquare")
S <- SRE(f,meuse,G)
g <- ggplot()+
    geom_point(data=data.frame(meuse),aes(x,y))
show_basis(g,S@basis) + coord_fixed()
