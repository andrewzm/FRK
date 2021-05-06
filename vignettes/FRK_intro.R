## ----setup, include=FALSE, cache=FALSE----------------------------------------
library(knitr)
# set global chunk options
# opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold')
# options(formatR.arrow=TRUE,width=90)
knitr::opts_chunk$set(dpi=100)

## ----eval=TRUE,message=FALSE,warning=FALSE------------------------------------
library(sp)        # for defining points/polygons
library(ggplot2)   # for plotting
library(dplyr)     # for easy data manipulation
library(FRK)       # for carrying out FRK

## ----eval=TRUE----------------------------------------------------------------
opts_FRK$set("progress",FALSE)  # no progress bars
opts_FRK$set("parallel",0L)     # no parallelisation

## -----------------------------------------------------------------------------
data(meuse)            # load meuse data
print(class(meuse))    # print class of meuse data

## -----------------------------------------------------------------------------
coordinates(meuse) = ~x+y     # change into an sp object

## ----message=FALSE------------------------------------------------------------
set.seed(1)
GridBAUs1 <- auto_BAUs(manifold = plane(),    # 2D plane
                     cellsize = c(100,100),   # BAU cellsize
                     type = "grid",           # grid (not hex)
                     data = meuse,            # data around which to create BAUs
                     convex=-0.05,            # border buffer factor
                     nonconvex_hull=FALSE)    # convex hull

## -----------------------------------------------------------------------------
GridBAUs1$fs <- 1   # fine-scale variation at BAU level

## ----echo=FALSE,fig.cap="(a) Locations of the \\tt{meuse} data. (b) BAUs for Fixed Rank Kriging with the \\tt{meuse} dataset.\\label{fig:meuse}",fig.subcap=c("",""),fig.width=5,fig.height=4,out.width="0.5\\linewidth",fig.align='center'----
plot(NULL,NULL,xlim = c(178605,181390),ylim=c(329714,333611),asp=1,xlab="Easting (m)",ylab="Northing (m)")
plot(meuse,add=T)
plot(NULL,NULL,xlim = c(178605,181390),ylim=c(329714,333611),asp=1,xlab="Easting (m)",ylab="Northing (m)")
plot(as(GridBAUs1,"SpatialPolygons"),add=T)

## ----message=FALSE------------------------------------------------------------
G <- auto_basis(manifold = plane(),   # 2D plane
                data = meuse,         # meuse data
                nres = 2,             # number of resolutions
                type = "Gaussian",    # type of basis function
                regular = 1)          # place regularly in domain

## ----message=FALSE, fig.cap="Basis functions automatically generated for the meuse dataset with 2 resolutions. The interpretation of the circles change with the domain and basis. For Gaussian functions on the plane, each circle is centred at the basis function centre, and has a radius equal to $1\\sigma$. Type \\tt{help(auto\\_basis)} for details.\\label{fig:basis}",fig.height=6,fig.width=6,out.width="0.5\\linewidth",fig.align='center',fig.pos="t"----
show_basis(G) +             # illustrate basis functions
    coord_fixed() +         # fix aspect ratio
    xlab("Easting (m)") +   # x-label
    ylab("Northing (m)")    # y-label

## -----------------------------------------------------------------------------
f <- log(zinc) ~ 1    # formula for SRE model

## ----results='hide', message=FALSE--------------------------------------------
S <- SRE(f = f,                  # formula
         data = list(meuse),     # list of datasets
         BAUs = GridBAUs1,       # BAUs
         basis = G,              # basis functions
         est_error = TRUE,       # estimation measurement error
         average_in_BAU = FALSE) # do not average data over BAUs

## ----message=FALSE,results='hide',cache=FALSE,fig.cap="Convergence of the EM algorithm when using \\tt{FRK} with the \\tt{meuse} dataset.\\label{fig:EM}",fig.height=4,fig.width=5,fig.align='center'----
S <- SRE.fit(S,                # SRE model
             n_EM = 10,        # max. no. of EM iterations
             tol = 0.01,       # tolerance at which EM is assumed to have converged
             print_lik=TRUE)   # print log-likelihood at each iteration

## -----------------------------------------------------------------------------
GridBAUs1 <- predict(S, obs_fs = FALSE)

## -----------------------------------------------------------------------------
BAUs_df <- as(GridBAUs1,"data.frame")

## -----------------------------------------------------------------------------
g1 <- ggplot() +                          # Use a plain theme
    geom_tile(data=BAUs_df ,                  # Draw BAUs
                 aes(x,y,fill=mu),      # Colour <-> Mean
                 colour="light grey") +          # Border is light grey
    scale_fill_distiller(palette="Spectral",     # Spectral palette
                         name="pred.") +         # legend name
    geom_point(data=data.frame(meuse),           # Plot data
               aes(x,y,fill=log(zinc)),          # Colour <-> log(zinc)
               colour="black",                   # point outer colour
               pch=21, size=3) +                 # size of point
    coord_fixed() +                              # fix aspect ratio
    xlab("Easting (m)") + ylab("Northing (m)") + # axes labels
    theme_bw()


g2 <- ggplot() +                          # Similar to above but with s.e.
    geom_tile(data=BAUs_df,
                 aes(x,y,fill=sqrt(var)),
                 colour="light grey") +
    scale_fill_distiller(palette="BrBG",
                         name = "s.e.",
                         guide = guide_legend(title="se")) +
    coord_fixed() +
    xlab("Easting (m)") + ylab("Northing (m)") + theme_bw()

## ----echo=FALSE,fig.cap="Inference at the BAU level using FRK with the \\tt{meuse} dataset. (a) FRK prediction. (b) FRK prediction standard error.\\label{fig:PredictionBAU}",fig.width=6,fig.height=7.5,out.width="0.5\\linewidth",fig.subcap=c('',''),fig.pos="t"----
plot(g1)
plot(g2)

## -----------------------------------------------------------------------------
Pred_regions <- auto_BAUs(manifold = plane(),      # model on the 2D plane
                          cellsize = c(600,600),   # choose a large grid size
                          type = "grid",           # use a grid (not hex)
                          data = meuse,            # the dataset on which to center cells
                          convex=-0.05,            # border buffer factor
                          nonconvex_hull=FALSE)    # convex hull

## -----------------------------------------------------------------------------
Pred_regions <- predict(S, newdata = Pred_regions)    # prediction polygons

## ----echo=FALSE,fig.cap="Prediction and prediction standard error obtained with FRK from the \\tt{meuse} dataset over arbitrary polygons. Both quantities are logs of ppm.\\label{fig:PredictionPolygon}",fig.subcap=c("",""),fig.width=6,fig.height=7.5,out.width="0.5\\linewidth",fig.pos="t"----
Pred_regions_df <- SpatialPolygonsDataFrame_to_df(
    sp_polys = as(Pred_regions, "SpatialPolygonsDataFrame"),
    vars = c("mu","var"))

g1 <- ggplot() +
    geom_polygon(data=Pred_regions_df,
                 aes(x,y,fill=mu,group=id),
                 colour="light grey") +
    scale_fill_distiller(palette="Spectral",name="pred.") +
    geom_point(data=data.frame(meuse),
               aes(x,y,fill=log(zinc)),
               colour="black",
               pch=21, size=3) +
    coord_fixed() +
    xlab("Easting (m)") + ylab("Northing (m)") + theme_bw()

g2 <- ggplot() +
    geom_polygon(data=Pred_regions_df,
                 aes(x,y,fill=sqrt(var),group=id),
                 colour="light grey") +
    scale_fill_distiller(palette="BrBG",
                         guide = guide_legend(title="s.e.")) +
    coord_fixed() +
    xlab("Easting (m)") + ylab("Northing (m)") + theme_bw()

plot(g1)
plot(g2)

## -----------------------------------------------------------------------------
data(meuse)
meuse[1:10,"zinc"] <- NA

## -----------------------------------------------------------------------------
meuse2 <- subset(meuse,!is.na(zinc))
meuse2 <- meuse2[,c("x","y","zinc")]
coordinates(meuse2) <- ~x+y

## -----------------------------------------------------------------------------
meuse$zinc <- NULL
coordinates(meuse) <- c("x","y")
meuse.grid2 <- BAUs_from_points(meuse)

## ----eval=TRUE----------------------------------------------------------------
data(AIRS_05_2003)                                          ## Load data

## -----------------------------------------------------------------------------
AIRS_05_2003 <-
    dplyr::filter(AIRS_05_2003,day %in% 1:3) %>%    # only first three days
    dplyr::mutate(std=co2std) %>%                   # change std to have suitable name
    dplyr::select(lon,lat,co2avgret,std)            # select columns we actually need
coordinates(AIRS_05_2003) = ~lon+lat                # change into an sp object
proj4string(AIRS_05_2003) =
            CRS("+proj=longlat +ellps=sphere")      # unprojected coordinates on sphere

## ----eval=TRUE, warning=FALSE,results='hide',message=FALSE--------------------
isea3h_sp_poldf <- auto_BAUs(manifold   = sphere(),  # model on sphere
                             isea3h_res = 6,         # isea3h resolution 6 BAUs
                             type = "hex",           # hexagonal grid
                             data = AIRS_05_2003)    # remove BAUs where there is not data
isea3h_sp_poldf$fs = 1                               # fine-scale component

## ----eval=TRUE, results='hide'------------------------------------------------
G <- auto_basis(manifold = sphere(), # basis functions on the sphere
                data=AIRS_05_2003,   # AIRS data
                nres = 2,            # number of resolutions
                type = "bisquare")   # bisquare function


## ----echo=FALSE,message = FALSE, fig.cap="BAUs and basis functions used in modelling and predicting with the \\tt{AIRS} data. (a) BAUs are the ISEA3H hexagons at resolution 5. (b) Basis function centroids constructed using the function \\tt{auto\\_basis}.\\label{fig:sphere_BAUs}",fig.subcap=c("",""),fig.width=6,fig.height=6,out.width="0.5\\linewidth",fig.pos="t!",message=FALSE,dev = 'png',dev = 'png'----
data("isea3h")
ggplot(subset(isea3h,res==5 & centroid==0)) +
          geom_path(aes(lon,lat,group=id)) +
          coord_map("ortho") +
          xlab("lon (deg)") +
          ylab("lat (deg)") + theme_bw()
show_basis(G,draw_world()) +
    coord_fixed(ratio = 2) +
          xlab("lon (deg)") +
          ylab("lat (deg)") + theme_bw()

## ----cache=FALSE,eval=TRUE,message=FALSE,results='hide'-----------------------
f <- co2avgret ~ lat + 1         # formula for fixed effects
S <- SRE(f = f,                  # formula for fixed effects
         list(AIRS_05_2003),     # list of data objects
         basis = G,              # basis functions
         BAUs = isea3h_sp_poldf, # BAUs
         est_error = FALSE,      # do not estimate meas. error
         average_in_BAU = TRUE)  # summarise data

S <- SRE.fit(S,                # SRE model
             n_EM = 1,         # max. no. of EM iterations
             tol = 0.01,       # tolerance at which EM is assumed to have converged
             print_lik=FALSE)  # do not print log-likelihood at each iteration

## ----eval=TRUE----------------------------------------------------------------
isea3h_sp_poldf <- predict(S)          # fs variation is in the observation model

## ----echo=FALSE,results='hide',message=FALSE----------------------------------
X <- SpatialPolygonsDataFrame_to_df(sp_polys = isea3h_sp_poldf,
                                    vars = c("mu","var"))
mumin <- quantile(X$mu,0.01)
mumax <- quantile(X$mu,0.99)

## ----echo=FALSE,fig.width=7,fig.height=3.5,out.width="\\linewidth",fig.pos="t!",fig.cap="CO$_2$ mole-fraction readings in ppm from the \\tt{AIRS}.\\label{fig:AIRSresults1}",fig.align="centre",dev = 'png',dpi=70----
g1 <- (EmptyTheme() +
           geom_point(data=data.frame(AIRS_05_2003),
                      aes(lon,lat,
                          colour=pmin(pmax(
                              co2avgret,mumin),
                              mumax)),
                      pch=46) +
           scale_colour_distiller(palette="Spectral",
                                  guide_legend(title="co2 (ppm)")
           ) +
           coord_map("mollweide")) %>%
    draw_world(inc_border=TRUE) +
          xlab("lon (deg)") +
          ylab("lat (deg)")
print(g1)

## ----echo=FALSE,fig.width=7,fig.height=3.5,out.width="\\linewidth",fig.pos="t!",fig.cap="Prediction of $\\Yvec_P$ in ppm following FRK on the \\tt{AIRS} data.\\label{fig:AIRSresults2}",fig.align="centre",dev = 'png',dpi=70----
g2 <- (EmptyTheme() +
           geom_polygon(data=X,
                        aes(lon,lat,fill=pmin(pmax(mu,mumin),mumax),group=id))+
           scale_fill_distiller(palette="Spectral",name="pred. (ppm)",limits=c(mumin,mumax)) +
           coord_map("mollweide")) %>%
    draw_world(inc_border=TRUE)+
          xlab("lon (deg)") +
          ylab("lat (deg)")
print(g2)

## ----echo=FALSE,fig.keep=TRUE,fig.width=7,fig.height=3.5,out.width="\\linewidth",fig.pos="t!",fig.cap="Prediction standard error of $\\Yvec_P$ in ppm following FRK on the \\tt{AIRS} data.\\label{fig:AIRSresults3}",fig.align="centre",dev = 'png',dpi=70----
X$se <- pmax(sqrt(X$var),0.32)
g3 <- (EmptyTheme() +
           geom_polygon(data=X,
                        aes(lon,lat,fill=se,group=id))+
           scale_fill_distiller(palette="BrBG",name="s.e. (ppm)") +
           coord_map("mollweide")) %>%
    draw_world(inc_border=TRUE)+
          xlab("lon (deg)") +
          ylab("lat (deg)")

print(g3)

## -----------------------------------------------------------------------------
library(spacetime)

## ----message=FALSE------------------------------------------------------------
data("NOAA_df_1990")             # load data
Tmax <- subset(NOAA_df_1990,     # subset the data
              month %in% 7 &     # May to July
              year == 1993)      # year of 1993

## ----echo=FALSE,fig.align="center",fig.height=4,fig.width=6,fig.cap="Station locations from which the maximum temperature readings in the \\tt{NOAA} dataset were obtained.\\label{fig:stat_locs}",out.width="0.6\\linewidth"----
(ggplot() + geom_point(data=Tmax,aes(lon,lat),size=0.5)) %>%
    draw_world() + coord_fixed(xlim = c(-115,-60), ylim = c(10,60)) + theme_bw()

## -----------------------------------------------------------------------------
Tmax <- within(Tmax,
               {time = as.Date(paste(year,month,day,sep="-"))})  # create Date field

## -----------------------------------------------------------------------------
STObj <- stConstruct(x = Tmax,                    # dataset
                 space = c("lon","lat"),          # spatial fields
                 time="time",                     # time field
                 interval=TRUE)                   # time reflects an interval

## -----------------------------------------------------------------------------
STObj$std <- 2

## ----warning=FALSE------------------------------------------------------------
grid_BAUs <- auto_BAUs(manifold=STplane(),    # spatio-temporal process on the plane
                       data=STObj,            # data
                       cellsize = c(1,1,1),   # BAU cell size
                       type="grid",           # grid or hex?
                       convex=-0.1,           # parameter for hull construction
                       tunit="days",          # time unit
                       nonconvex_hull=FALSE)  # convex hull
grid_BAUs$fs = 1                       # fine-scale variation

## -----------------------------------------------------------------------------
G_spatial <- auto_basis(manifold = plane(),          # spatial functions on the plane
                        data=as(STObj,"Spatial"),    # remove the temporal dimension
                        nres = 1,                    # three resolutions
                        type = "bisquare",           # bisquare basis functions
                        regular = 1)                 # regular basis functions

## ----warning=FALSE------------------------------------------------------------
print(head(grid_BAUs@time))                                # show time indices
G_temporal <- local_basis(manifold = real_line(),          # functions on the real line
                           type = "Gaussian",              # Gaussian functions
                           loc = matrix(seq(2,28,by=4)),   # locations of functions
                           scale = rep(3,7))               # scales of functions

## ----message=FALSE------------------------------------------------------------
basis_s_plot <- show_basis(G_spatial) + xlab("lon (deg)") + ylab("lat (deg)")
basis_t_plot <- show_basis(G_temporal) + xlab("time index") + ylab(expression(phi(t)))

## ----echo=FALSE,fig.height=4,fig.width=4,fig.subcap=c("",""),fig.cap="Spatial and temporal basis functions used to construct the spatio-temporal basis functions. (a)  Spatial support of the bisquare spatial basis functions. (b) The temporal basis functions.\\label{fig:STbasis}",out.width="0.5\\linewidth"----
print(basis_s_plot)
print(basis_t_plot)

## -----------------------------------------------------------------------------
G <- TensorP(G_spatial,G_temporal)         # take the tensor product

## ----SRE,results='hide',cache=FALSE-------------------------------------------
 f <- z ~ 1 + lat                 # fixed effects part
 S <- SRE(f = f,                  # formula
          data = list(STObj),     # data (can have a list of data)
          basis = G,              # basis functions
          BAUs = grid_BAUs,       # BAUs
          est_error = FALSE)      # do not estimate measurement-error variance

 S <- SRE.fit(S,                # estimate parameters in the SRE model S
             n_EM = 1,          # maximum no. of EM iterations
             tol = 0.1,         # tolerance on log-likelihood
             print_lik=FALSE)   # print log-likelihood trace

 grid_BAUs <- predict(S, obs_fs = FALSE)

## ----message=FALSE,results='hide'---------------------------------------------
 analyse_days <- c(1,4,8,12,16,20)  # analyse only a few days
 df_st <- lapply(analyse_days,      # for each day
        function(i)
            as(grid_BAUs[,i],"data.frame") %>%
            cbind(day = i))         # add day number to df
 df_st <- do.call("rbind",df_st)    # append all dfs together

## ----echo=FALSE,fig.cap="Spatio-temporal FRK prediction of \\tt{Tmax} on the plane in degrees Fahrenheit within a domain enclosing the region of interest for six selected days spanning the temporal window of the data: 01 July 1993 -- 20 July 2003.\\label{fig:FRK_pred1}",fig.height=8,fig.width=16,fig.pos="t!"----
ggplot() +                          # Similar to above but with s.e.
    geom_tile(data=df_st,
                 aes(lon,lat,fill=mu),
                 colour="light grey") +
      geom_point(data=filter(Tmax,day %in% c(1,4,8,12,16,20)),           # Plot data
               aes(lon,lat,fill=z),          # Colour <-> log(zinc)
               colour="black",                   # point outer colour
               pch=21, size=3) +                 # size of point
    scale_fill_distiller(palette="Spectral",
                         guide = guide_legend(title="pred. (degF)")) +
    coord_fixed() +
    xlab("lon (deg)") + ylab("lat (deg)") +
    facet_wrap(~day) + theme_bw()

## ----echo=FALSE,fig.cap="Spatio-temporal FRK prediction standard  error of  \\tt{Tmax} on the plane in degrees Fahrenheit within a domain enclosing the region of interest for the same six days selected in Fig.~\\ref{fig:FRK_pred1} and spanning the temporal window of the data, 01 July 1993 -- 20 July 2003.\\label{fig:FRK_pred2}",fig.height=8,fig.width=16,fig.pos="t!"----
ggplot() +                          # Similar to above but with s.e.
    geom_tile(data=df_st,
                 aes(lon,lat,fill=sqrt(var)),
                 colour="light grey") +
      scale_fill_distiller(palette="BrBG",
                         trans="reverse",
                         guide = guide_legend(title="s.e. (degF)")) +
    coord_fixed() +
    xlab("lon (deg)") + ylab("lat (deg)") +
    facet_wrap(~day) + theme_bw()

## -----------------------------------------------------------------------------
proj4string(STObj) <- "+proj=longlat +ellps=sphere"

## ----FRK2,cache=FALSE---------------------------------------------------------
grid_BAUs <- auto_BAUs(manifold=STsphere(),       # spatio-temporal process on the sphere
                       data=STObj,                # data
                       cellsize = c(1,1,1),       # BAU cell size
                       type="grid",               # grid or hex?
                       convex=-0.1,               # parameter for hull construction
                       tunit="days")              # time unit

 G_spatial <- auto_basis(manifold = sphere(),      # spatial functions on the plane
                        data=as(STObj,"Spatial"),  # remove the temporal dimension
                        nres = 2,                  # two resolutions of DGG
                        type = "bisquare",         # bisquare basis functions
                        prune=15,                  # prune basis functions
                        isea3h_lo = 4)             # but remove those lower than res 4

## ----message=FALSE,echo=FALSE,fig.cap="Basis functions for FRK on the sphere with the \\tt{NOAA} dataset using two ISEA3H DGGs for location parameters of the basis functions.\\label{fig:basis_USA}",fig.height=9,fig.width=9,fig.pos="t!",out.width="0.7\\linewidth",fig.align="center"----
draw_world(show_basis(G_spatial,ggplot())) + coord_map("ortho",orientation = c(35,-100,0)) + xlab("lon (deg)") + ylab("lat (deg)") + theme_bw()

## ----eval=TRUE----------------------------------------------------------------
data(AIRS_05_2003)   # load AIRS data

## ----eval=TRUE----------------------------------------------------------------
set.seed(1)
AIRS_05_2003 <- mutate(AIRS_05_2003,           # take the data
                       std=co2std) %>%         # rename std
                sample_n(20000)                # sample 20000 points

## -----------------------------------------------------------------------------
AIRS_05_2003 <- within(AIRS_05_2003,
               {time = as.Date(paste(year,month,day,sep="-"))})  # create Date field

## -----------------------------------------------------------------------------
STObj <- stConstruct(x = AIRS_05_2003,            # dataset
                 space = c("lon","lat"),          # spatial fields
                 time ="time",                    # time field
                 crs = CRS("+proj=longlat +ellps=sphere"),  # CRS
                 interval=TRUE)                   # time reflects an interval

## ----eval=TRUE,cache=FALSE----------------------------------------------------
## Prediction (BAU) grid
grid_BAUs <- auto_BAUs(manifold=STsphere(),         # space-time field on sphere
                             data=time(STObj),      # temporal part of the data
                             cellsize = c(5,5,1),   # cellsize (5 deg x 5 deg x 1 day)
                             type="grid",           # grid (not hex)
                             tunit = "days")        # time spacing in days
grid_BAUs$fs = 1

## ----echo=FALSE,fig.subcap=c("",""),fig.cap="Gridded BAUs on the sphere used for modelling and predicting with the \\tt{AIRS} data. (a) BAUs constructed when supplying only the temporal indices of the data (the entire sphere is covered with BAUs within the specified time period). (b) BAUs constructed when supplying the entire dataset. The view of the sphere is from the bottom; in this case there is no data below $60^{\\circ}$S and thus BAUs have been omitted from this region.\\label{fig:sphere_grid_BAUs}",fig.width=4,fig.height=4,out.width="0.5\\linewidth",fig.pos="t!",message=FALSE----

grid_BAUs2 <- auto_BAUs(manifold=STsphere(),  # space-time field on sphere
                             data=STObj,            # data
                             cellsize = c(5,5,1),   # cellsize (5 deg x 5 deg x 1 day)
                             type="grid",           # grid (not hex)
                             tunit = "days")        # time spacing in days

X <- SpatialPolygonsDataFrame_to_df(grid_BAUs[,1],"n")
X2 <- SpatialPolygonsDataFrame_to_df(grid_BAUs2[,1],"n")
ggplot(X) +
    geom_polygon(aes(lon,lat,group=id),colour="black",fill="white") +
    coord_map("ortho",orientation = c(-145,125,25)) +
    xlab("lon (deg)") + ylab("lat (deg)") + theme_bw()
ggplot(X2) +
    geom_polygon(aes(lon,lat,group=id),colour="black",fill="white") +
    coord_map("ortho",orientation = c(-145,125,25)) +
    xlab("lon (deg)") + ylab("lat (deg)") + theme_bw()

## ----eval=TRUE----------------------------------------------------------------
G_spatial <- auto_basis(manifold = sphere(),      # functions on sphere
                        data=as(STObj,"Spatial"), # collapse time out
                        nres = 1,                 # use three DGGRID resolutions
                        prune= 15,                 # prune basis functions
                        type = "bisquare",        # bisquare basis functions
                        subsamp = 2000,          # use only 2000 data points for pruning
                        isea3h_lo = 2)            # start from isea3h res 2

G_temporal <- local_basis(manifold=real_line(),      # functions on real line
                          loc = matrix(c(2,7,12)),   # location parameter
                          scale = rep(3,3),          # scale parameter
                          type = "Gaussian")
G_spacetime <- TensorP(G_spatial,G_temporal)

## ----eval=TRUE,message=FALSE,results='hide',cache=FALSE-----------------------
f <- co2avgret ~ lat +1           # formula for fixed effects
S <- SRE(f = f,                   # formula
         data = list(STObj),      # spatio-temporal object
         basis = G_spacetime,     # space-time basis functions
         BAUs = grid_BAUs,        # space-time BAUs
         est_error = FALSE,       # do not estimate measurement error
         average_in_BAU = TRUE)   # average data that fall inside BAUs

S <- SRE.fit(S,                   # SRE model
             n_EM = 1,            # max. EM iterations
             tol = 0.01)          # convergence criteria

grid_BAUs <- predict(S, obs_fs = TRUE,         # fs variation is in obs. model
                     pred_time = c(4L,8L,12L)) # predict only at select days

## ----echo=FALSE,message=FALSE,results='hide'----------------------------------
X <- lapply(1:length(time(grid_BAUs)),
            function(i) {
                SpatialPolygonsDataFrame_to_df(sp_polys = grid_BAUs[,i],
                                    vars = c("mu","var")) %>%
            mutate(t = as.vector(grid_BAUs@time[i]))})
X <- do.call("rbind",X)
mumin <- min(X$mu)
mumax <- max(X$mu)

## ----echo=FALSE, fig.keep=TRUE,fig.cap="CO$_2$ readings taken from the \\tt{AIRS} on the 04, 08 and 12 May 2003 in ppm. \\label{fig:FRK_AIRS_ST1}",fig.align="center",out.width="\\linewidth",fig.height=3,fig.width=16,fig.pos="t!",dev = 'png',dpi=70----
g1 <- (EmptyTheme() +
           geom_point(data=dplyr::filter(AIRS_05_2003,day %in% c(4,8,12)),
                      aes(lon,lat,
                          colour=pmin(pmax(
                              co2avgret,mumin),
                              mumax)),
                      size=0.5) +
           facet_grid(~day)+
           scale_colour_distiller(palette="Spectral",
                                  guide_legend(title="co2 (ppm)")) +
           coord_map("mollweide") +
            xlab("lon (deg)") +
            ylab("lat (deg)")) %>%
    draw_world(inc_border=TRUE)
print(g1)

## ----echo=FALSE, fig.keep=TRUE,fig.cap="Prediction of $\\Yvec_P$ in ppm on 04, 08, and 12 May 2003 obtained with \\pkg{FRK} on the \\tt{AIRS} data. \\label{fig:FRK_AIRS_ST2}",fig.align="center",out.width="\\linewidth",fig.height=3,fig.width=16,fig.pos="t!",dev = 'png',dpi=70----

g2 <-  (EmptyTheme() +
            geom_polygon(data=filter(X,(t %in% c(4,8,12)) & abs(lon) < 175),
                         aes(lon,lat,fill=mu,group=id))+
           scale_fill_distiller(palette="Spectral",name="pred. (ppm)") +
           coord_map("mollweide") +
           facet_grid(~t) +
            xlab("lon (deg)") +
            ylab("lat (deg)")) %>%
    draw_world(inc_border=FALSE)

print(g2)

## ----echo=FALSE, fig.keep=TRUE,fig.cap="Prediction standard error of $\\Yvec_P$ in ppm on 04, 08 and 12 May 2003 obtained with \\pkg{FRK} on the \\tt{AIRS} data. \\label{fig:FRK_AIRS_ST3}",fig.align="center",out.width="\\linewidth",fig.height=3,fig.width=16,fig.pos="t!",dev = 'png',dpi=70----

X$se <- sqrt(X$var)
g3 <-  (EmptyTheme() +
            geom_polygon(data=filter(X,(t %in% c(4,8,12)) & abs(lon) < 175),
                         aes(lon,lat,fill=se,group=id))+
           scale_fill_distiller(palette="Spectral",name="s.e. (ppm)") +
           coord_map("mollweide") +
           facet_grid(~t) +
            xlab("lon (deg)") +
            ylab("lat (deg)")) %>%
    draw_world(inc_border=FALSE)

print(g3)

## ----echo=FALSE,include=FALSE,warning=FALSE-----------------------------------
# Generate observations with large spatial support
data(meuse.grid)
data(meuse)

meuse_pols <- NULL
offset <- 150
for(i in 1:nrow(meuse)) {
    this_meuse <- meuse[i,]
    meuse_pols <- rbind(meuse_pols,
                        data.frame(x = c(this_meuse$x - offset,
                                         this_meuse$x + offset,
                                         this_meuse$x + offset,
                                         this_meuse$x - offset),
                                   y = c(this_meuse$y - offset,
                                         this_meuse$y - offset,
                                         this_meuse$y + offset,
                                         this_meuse$y + offset),
                                   id = i,
                                   zinc = this_meuse$zinc))
}
meuse_pols <- df_to_SpatialPolygons(meuse_pols,coords=c("x","y"),keys="id",proj = CRS())
meuse_pols <- SpatialPolygonsDataFrame(meuse_pols,data.frame(row.names = row.names(meuse_pols),zinc=meuse$zinc))
coordnames(meuse_pols) <- c("x","y")
coordinates(meuse) = ~x + y
meuse_pols$zinc <- exp(over(meuse_pols,GridBAUs1)$mu)

## ----message=FALSE,echo=FALSE,cache=FALSE,results='hide',warning=FALSE--------
set.seed(1)
GridBAUs2 <- auto_BAUs(manifold = plane(),     # 2D plane
                     cellsize = c(100,100),   # BAU cellsize
                     type = "grid",           # grid (not hex)
                     data = meuse,            # data around which to create BAUs
                     convex=-0.05,            # border buffer factor
                     nonconvex_hull = 0)               # convex hull
GridBAUs2$fs <- 1   # fine-scale variation at BAU level
G <- auto_basis(manifold = plane(),   # 2D plane
                data=meuse,           # meuse data
                nres = 2,             # number of resolutions
                type = "Gaussian",    # type of basis function
                regular = 1,prune=1)          # place regularly in domain
f <- log(zinc) ~ 1    # formula for SRE model
meuse_pols$std <- 1
S <- SRE(f = f,                # formula
         data = list(meuse_pols),   # list of datasets
         BAUs = GridBAUs2,      # BAUs
         basis = G,            # basis functions
         est_error=TRUE)       # estimation measurement error
S <- SRE.fit(S,                # SRE model
             n_EM = 4,         # max. no. of EM iterations
             tol = 0.01,       # tolerance at which EM is assumed to have converged
             print_lik=FALSE)   # print log-likelihood at each iteration
GridBAUs2 <- predict(S, obs_fs = FALSE)
BAUs_df <- as(GridBAUs2,"data.frame")
Obs_df <- SpatialPolygonsDataFrame_to_df(sp_polys = meuse_pols,   # BAUs to convert
                                          vars = c("zinc"))  # fields to extract
g1 <- ggplot() +                          # Use a plain theme
    geom_tile(data=BAUs_df ,                  # Draw BAUs
                 aes(x,y,fill=mu),      # Colour <-> Mean
                 colour="light grey") +          # Border is light grey
    scale_fill_distiller(palette="Spectral",name="pred.")  +  # Spectral palette
    coord_fixed() +                              # fix aspect ratio
    xlab("Easting (m)") + ylab("Northing (m)") + # axes labels
    theme_bw()

g2 <- ggplot() +                          # Similar to above but with s.e.
    geom_tile(data=BAUs_df,
                 aes(x,y,fill=sqrt(var)),
                 colour="light grey") +
    scale_fill_distiller(palette="BrBG",
                         guide = guide_legend(title="s.e.")) +
    coord_fixed() +
    xlab("Easting (m)") + ylab("Northing (m)") +
    geom_path(data=Obs_df,
                 aes(x,y,group=id),
                 colour="black",
                 alpha = 0.5) + theme_bw()


## ----echo=FALSE,fig.cap="Prediction and prediction standard error obtained with FRK using the \\tt{meuse} dataset where each observation is assuming to have a spatial footprint of 300 m $\\times$ 300m. (a) FRK prediction at the BAU level. (b) FRK prediction standard error at the BAU level. The black hexagons outline the spatial footprints of the data.\\label{fig:meuse_large}",fig.width=6,fig.height=7.5,out.width="0.5\\linewidth",fig.subcap=c("",""),fig.pos="t",dpi=10----
plot(g1)
plot(g2)

## -----------------------------------------------------------------------------
set.seed(1)
N <- 50
sim_process <- expand.grid(x = seq(0.005,0.995,by=0.01),       # x grid
                           y = seq(0.001,0.995,by=0.01)) %>%   # y grid
    mutate(proc = cos(x*40)*cos(y*3) + 0.3*rnorm(length(x)))   # anisotropic function

sim_data <- sample_n(sim_process,1000) %>%                     # sample data from field
    mutate(z = proc + 0.1*rnorm(length(x)),                    # add noise
           std = 0.1,                                          # with 0.1 std
           x = x + runif(1000)*0.001,                          # jitter x locations
           y = y + runif(1000)*0.001)                          # jitter y locations
coordinates(sim_data) = ~x + y                                 # change into SpatialPoints

## ----echo=FALSE,eval=TRUE,fig.cap="FRK with anisotropic fields. (a) Simulated process. (b) Observed data. \\label{fig:aniso1}",fig.width=6,fig.height=5,out.width="0.5\\linewidth",fig.subcap=c('',''),fig.pos="t"----
 g1 <-ggplot() +
     scale_fill_distiller(palette="Spectral",name="Y") +
     geom_tile(data=sim_process,aes(x,y,fill=proc))+
     coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
      xlab(expression(s[1])) + ylab(expression(s[2])) + theme_bw()

g2 <- ggplot() +
    scale_fill_distiller(palette="Spectral") +
    geom_point(data=data.frame(sim_data),
               aes(x,y,fill=z),
               colour="black",
               pch=21, size=2) +
    coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
      xlab(expression(s[1])) + ylab(expression(s[2])) + theme_bw()

  print(g1)
  print(g2)

## ----eval=TRUE----------------------------------------------------------------
scaler <- diag(c(4,1))                                    # scale x by 4
asymm_measure <- new("measure",                           # new measure object
                      dist=function(x1,x2=x1)                # new distance function
                            FRK:::distR(x1 %*% scaler,    # scaling of first point
                                        x2 %*% scaler),   # scaling of second point
                      dim=2L)                             # in 2D

## -----------------------------------------------------------------------------
TwoD_manifold <- plane()                 # Create R2 plane
TwoD_manifold@measure <- asymm_measure   # Assign measure

## -----------------------------------------------------------------------------
basis_locs <- seq(0,1,length=14) %>%                 # x locations
   expand.grid(seq(0,1,length=5)) %>%               # y locations
   as.matrix()                                      # convert to matrix
G <-  local_basis(manifold = TwoD_manifold,          # 2D plane
                 loc=basis_locs,                    # basis locations
                 scale=rep(0.4,nrow(basis_locs)),   # scale parameters
                 type="bisquare")                   # type of function

## ----echo=FALSE,eval=TRUE,fig.cap="Basis function 23 of the 75 constructed to fit an anisotropic spatial field. Anisotropy is obtained by changing the \\tt{measure} object of the manifold on which the basis function is constructed.\\label{fig:anisobasis}",fig.width=6,fig.height=6,out.width="0.5\\linewidth",fig.pos="t",fig.align="center"----
S <- eval_basis(G,as.matrix(sim_process[c("x","y")]))
sim_process$S <- S[,23]
ggplot() +
    geom_tile(data=sim_process,aes(x,y,fill=S)) +
    coord_fixed() + scale_fill_distiller(palette = "Spectral",guide =
                                            guide_legend(title=expression(phi[23](s)))) +
       xlab(expression(s[1])) + ylab(expression(s[2])) + theme_bw()

## ----echo=FALSE,cache=FALSE,message=FALSE,results='hide'----------------------
 ## Prediction (BAU) grid
 grid_BAUs <- auto_BAUs(manifold=plane(),
                        data=sim_data,
                        cellsize = c(0.02,0.02),
                        type="grid",
                        convex = -0.1,
                        nonconvex_hull=FALSE)
 grid_BAUs$fs = 1

  f <- z ~ 1
 S <- SRE(f = f,
          data = list(sim_data),
          basis = G,
          BAUs = grid_BAUs,
          est_error = FALSE,
          average_in_BAU = FALSE)

  S <- SRE.fit(S,
              n_EM = 4,
              tol = 0.01)

  grid_BAUs <- predict(S, obs_fs = TRUE)

   X <- as(grid_BAUs,"data.frame") %>%
     filter(x < 1.1 & x > -0.1 & y > -0.5 & y < 10.5)

 X$se <- sqrt(X$var)

g1 <-ggplot() +
    scale_fill_distiller(palette="Spectral",name = "pred.") +
    geom_tile(data=X,aes(x,y,fill=mu))+
    coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
      xlab(expression(s[1])) + ylab(expression(s[2]))

g2 <-ggplot() +
    scale_fill_distiller(palette="BrBG",name ="s.e.") +
    geom_tile(data=X,aes(x,y,fill=se))+
    coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
      xlab(expression(s[1])) + ylab(expression(s[2]))


## ----echo=FALSE,fig.subcap=c("",""),fig.cap="FRK using data generated by an anisotropic field. (a) FRK prediction. (b) FRK prediction standard error.\\label{fig:aniso2}",fig.width=6,fig.height=5,out.width="0.5\\linewidth",fig.pos="t"----
print(g1)
print(g2)

