# SpatialPolygonsDataFrame to df

Convert `SpatialPolygonsDataFrame` or `SpatialPixelsDataFrame` object to
data frame.

## Usage

``` r
SpatialPolygonsDataFrame_to_df(sp_polys, vars = names(sp_polys))
```

## Arguments

- sp_polys:

  object of class `SpatialPolygonsDataFrame` or `SpatialPixelsDataFrame`

- vars:

  variables to put into data frame (by default all of them)

## Details

This function is mainly used for plotting `SpatialPolygonsDataFrame`
objects with `ggplot` rather than `spplot`. The coordinates of each
polygon are extracted and concatenated into one long data frame. The
attributes of each polygon are then attached to this data frame as
variables that vary by polygon `id` (the rownames of the object).

## Examples

``` r
library(sp)
library(ggplot2)
opts_FRK$set("parallel",0L)
df <- data.frame(id = c(rep(1,4),rep(2,4)),
                 x = c(0,1,0,0,2,3,2,2),
                 y=c(0,0,1,0,0,1,1,0))
pols <- df_to_SpatialPolygons(df,"id",c("x","y"),CRS())
polsdf <- SpatialPolygonsDataFrame(pols,data.frame(p = c(1,2),row.names=row.names(pols)))
df2 <- SpatialPolygonsDataFrame_to_df(polsdf)
if (FALSE) ggplot(df2,aes(x=x,y=y,group=id)) + geom_polygon() # \dontrun{}
```
