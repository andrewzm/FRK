library(ggplot2)
library(doMC)
devtools::load_all("~/Wollongong/pkgs/FRK",
                   export_all = FALSE)

## Data obtained from http://webpages.sou.edu/~sahrk/dgg/isea.old/gen/isea3h.html

isea3h <- NULL
registerDoMC(7)
isea3h <- foreach(j = 0:9,.combine = "rbind") %dopar% {
    print(paste0("Processing resolution ",j))
    dggrid_gen_to_df(paste0("~/Desktop/isea3h",j,".gen"))
}
save(isea3h,file="./inst/extdata/isea3h.rda")
#ggplot(subset(isea3h,res==5)) + geom_point(aes(lon,lat,colour=centroid)) + coord_map("ortho")



