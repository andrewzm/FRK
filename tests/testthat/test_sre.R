library(sp)
data(meuse)
meuse$std <- sqrt(0.05066)
coordinates(meuse) = ~x+y # change into an sp object
f <- log(zinc) ~ 1

S <- SRE(f,meuse,radial_basis(plane(),loc = matrix(c(179024,329733),nrow=1),scale=600))
g <- ggplot()+
    geom_point(data=data.frame(meuse),aes(x,y))
show_basis(g,S@basis) + coord_fixed()
