#####################################################################
## Point-level data and predictions with the meuse dataset
#####################################################################
library("FRK")
library("sp")
library("ggplot2")
library("gridExtra")

## Initialise seed to ensure reproducibility
set.seed(1)

## Load meuse data and assume first 10 observations are missing
data("meuse")
meuse[1:10,"zinc"] <- NA

## Create a complete-data data frame
meuse2 <- subset(meuse,!is.na(zinc))
meuse2 <- meuse2[,c("x","y","zinc")]
coordinates(meuse2) <- ~x+y

## Create BAUs around all prediction and observation points (after removing data from field)
meuse$zinc <- NULL
coordinates(meuse) <- c("x","y")
meuse.grid2 <- BAUs_from_points(meuse)

## Run the function FRK with the data frame and created BAUs
f <- log(zinc) ~ 1 + sqrt(dist)
S <- FRK(f = f,
         data = list(meuse2),
         BAUs = meuse.grid2,
         regular = 0)
Pred <- predict(S, obs_fs = FALSE)

## Plot the predictions
BAUs_df <- as.data.frame(Pred)  # convert the BAUs to data frame
gFRK_pred <- ggplot(BAUs_df) + geom_point(aes(x,y,colour = mu)) +
  scale_colour_distiller(palette = "Spectral") +
  xlab("Easting (m)") + ylab("Northing (m)") + theme_bw()

## Plot the prediction standard errors: Red circles denote observed data
gFRK_se <-ggplot(BAUs_df) + geom_point(data = as.data.frame(meuse2),
                                       aes(x,y), col = "red", pch = 1, size = 2) +
  geom_point(aes(x,y,colour = sd)) +
  scale_colour_distiller(palette = "BrBG") +
  xlab("Easting (m)") + ylab("Northing (m)") + theme_bw()

grid.arrange(gFRK_pred,gFRK_se,nrow=1)
