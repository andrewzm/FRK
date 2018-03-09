#####################################################################
## BAU-level data and predictions with the meuse dataset
#####################################################################
library("FRK")
library("sp")
library("ggplot2")
library("gridExtra")

## Initialise seed to ensure reproducibility
set.seed(1)

## load the meuse data and convert to sp object
data("meuse")
coordinates(meuse) = ~ x + y

## load meuse.grid and convert to SpatialPixelsDataFrame
data("meuse.grid")
coordinates(meuse.grid) = ~ x + y
gridded(meuse.grid) = TRUE

## Remove any common fields in the data and the BAUs (in FRK all covariates
## need to be with the BAUs)
meuse$soil <- meuse$dist <- meuse$ffreq <- NULL

## formula for SRE
f <- log(zinc) ~ 1 + sqrt(dist)

## Run FRK
S <- FRK(f = f,                # formula
        data = list(meuse),    # data (just meuse)
        BAUs = meuse.grid,     # the BAUs
        regular = 0)           # irregularly placed basis functions

## Predict at all the BAUs
Pred <- predict(S, obs_fs = FALSE)  # Case 2 in paper (default)

## Plot results
BAUs_df <- as.data.frame(Pred)  # convert the BAUs to data frame
gFRK_pred <- ggplot(BAUs_df) + geom_tile(aes(x,y,fill = mu)) +
             scale_fill_distiller(palette = "Spectral") +
            xlab("Easting (m)") + ylab("Northing (m)") + theme_bw()
gFRK_se <-ggplot(BAUs_df) + geom_tile(aes(x,y,fill = sd)) +
            scale_fill_distiller(palette = "BrBG") +
            xlab("Easting (m)") + ylab("Northing (m)") + theme_bw()
grid.arrange(gFRK_pred,gFRK_se,nrow=1)


