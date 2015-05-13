library(dplyr)
library(ggplot2)
library(RCurl)
library(doMC)
library(zoo)

isea3h <- NULL
registerDoMC(7)
isea3h <- foreach(j = 0:8,.combine = "rbind") %dopar% {
    print(paste0("Processing resolution ",j))
    ## Data obtained from http://webpages.sou.edu/~sahrk/dgg/isea.old/gen/isea3h.html
    X <- read.table(paste0("~/Desktop/isea3h",j,".gen"),
                    sep=" ",fill=T,
                    header=F,
                    col.names = c("id","lon","lat")) %>%
        filter(!(id == "END")) %>%
        mutate(res = j,
               id = as.numeric(as.character(id)),
               centroid = as.numeric(!is.na(id))) %>%
        transform(id = na.locf(id))
    X
}
save(isea3h,file="./inst/extdata/isea3h.rda")
#ggplot(subset(isea3h,res==5)) + geom_point(aes(lon,lat,colour=centroid)) + coord_map("ortho")



