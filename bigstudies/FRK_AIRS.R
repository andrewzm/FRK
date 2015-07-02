### Load packages, including FRK
library(sp)
library(ggplot2)
library(dplyr)
library(FRK)

## Do not monitor progress and parallelise using 4 cores
opts_FRK$set("progress",TRUE)
opts_FRK$set("parallel",4L)

## Load data
data(AIRS_05_2003) 

## Use only first day, set std=co2std as field name and select the important columns
AIRS_05_2003 <- filter(AIRS_05_2003,day %in% 1:1) %>%
    mutate(std=co2std) %>%
    select(lon,lat,co2avgret,std)

## Convert to an "sp" object by assigning coordinates
coordinates(AIRS_05_2003) = ~lon+lat # change into an sp object

## Set projection
proj4string(AIRS_05_2003)=CRS("+proj=longlat")

## Set the Basic Aerial Units, these are the smallest units on which we will carry out inference
isea3h_sp_poldf <- auto_BAUs(manifold=sphere(), # Use the sphere manifold
                             data=NULL,         # Do not tune with data
                             cellsize = c(2,2), # Set grid to 5 lat/lons
                             type="grid")       # Arrange as grid (can be "hex")
isea3h_sp_poldf$fs = 1                          # Assign unit fine-scale variation to BAUs

## Set up basis functions for the model
G <- auto_basis(m = sphere(),               # Use sphere manifold
                data=AIRS_05_2003,          # Use AIRS data
                nres = 3,                   # Only use 3 resolutions
                prune=15,                   # Remove basis functions when sum exceeds 15
                type = "bisquare",          # Use bisquare basis (can be set to Gaussian)
                subsamp = 20000)            # Use only 20,000 data points for filtering out basis functions

## To plot the centres and radii of basis functions do
show_basis(draw_world(inc_border=TRUE),G) + coord_map("mollweide")

## Model we are going to fit
f <- co2avgret ~ lat + 1           # Assume fixed effect latitude + intercept
S <- SRE(f,                        # Construct the spatial random effects model
         list(AIRS_05_2003),       # Use only one dataset (can use multiple)
         G,                        # Basis functions to use
         isea3h_sp_poldf,          # BAUs for prediction
         est_error = FALSE)        # estimate observation error?

## Fit the model using the EM algorithm
S <- SRE.fit(S,                # Model object
             n_EM = 100,         # Number of EM iterations
             tol = 1e-5,       # Convergence tolerance
             print_lik=TRUE)   # Print likelihood ?


## Predict at BAU locations
isea3h_sp_poldf <- SRE.predict(S,
                               pred_locs = isea3h_sp_poldf,
                               use_centroid = TRUE) # Use centroid only for computation (quicker)


## Convert polygons to data frame for plotting
X <- SpatialPolygonsDataFrame_to_df(sp_polys = isea3h_sp_poldf,
                                    vars = c("mu","var"))

## Plot
g1 <- (EmptyTheme() +
           geom_polygon(data=X,
                        aes(lon,lat,fill=mu,group=id))+
           scale_fill_distiller(palette="Spectral",trans="reverse") +
           coord_map("mollweide")) %>%
    draw_world(inc_border=FALSE)

mumin <- min(X$mu)
mumax <- max(X$mu)

g2 <- (EmptyTheme() +
           geom_point(data=data.frame(AIRS_05_2003),
                      aes(lon,lat,
                          colour=pmin(pmax(
                              co2avgret,mumin),
                              mumax)),
                      size=2) +
           scale_colour_distiller(palette="Spectral",
                                  trans="reverse",
                                  guide_legend(title="co2")
           ) +
           coord_map("mollweide")) %>%
    draw_world(inc_border=TRUE)

print(g1)
#print(g2)

