library(dplyr)
library(ggplot2)
library(RCurl)

isea3h <- NULL

for(j in 0:5) {
    print(paste0("Processing resolution ",j))
    ## Data obtained from http://webpages.sou.edu/~sahrk/dgg/isea.old/gen/isea3h.html
    X <- read.table(paste0("~/Desktop/isea3h",j,".gen"),sep=" ",fill=T,header=F,col.names = c("id","lon","lat")) %>%
         filter(!(id == "END"))

    X$res = j
    X$centroid = 0
    for(i in 1:nrow(X)) {
        if(!(X[i,1] == "")) {
            current_id <- X[i,1]
            X$centroid[i] = 1
        }
        X[i,1] <- current_di
    }
    X$id <- as.numeric(X$id)
    isea3h <- rbind(isea3h,X)
}
save(isea3h,file="./inst/extdata/isea3h.rda")
#ggplot(subset(isea3h,res==5)) + geom_point(aes(lon,lat,colour=centroid)) + coord_map("ortho")



