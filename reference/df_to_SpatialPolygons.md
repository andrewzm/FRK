# Convert data frame to SpatialPolygons

Convert data frame to SpatialPolygons object.

## Usage

``` r
df_to_SpatialPolygons(df, keys, coords, proj)
```

## Arguments

- df:

  data frame containing polygon information, see details

- keys:

  vector of variable names used to group rows belonging to the same
  polygon

- coords:

  vector of variable names identifying the coordinate columns

- proj:

  the projection of the `SpatialPolygons` object. Needs to be of class
  `CRS`

## Details

Each row in the data frame `df` contains both coordinates and labels (or
keys) that identify to which polygon the coordinates belong. This
function groups the data frame according to `keys` and forms a
`SpatialPolygons` object from the coordinates in each group. It is
important that all rings are closed, that is, that the last row of each
group is identical to the first row. Since `keys` can be of length
greater than one, we identify each polygon with a new key by forming an
MD5 hash made out of the respective `keys` variables that in themselves
are unique (and therefore the hashed key is also unique). For lon-lat
coordinates use `proj = CRS("+proj=longlat +ellps=sphere")`.

## Examples

``` r
library(sp)
df <- data.frame(id = c(rep(1,4),rep(2,4)),
                 x = c(0,1,0,0,2,3,2,2),
                 y=c(0,0,1,0,0,1,1,0))
pols <- df_to_SpatialPolygons(df,"id",c("x","y"),CRS())
if (FALSE) plot(pols) # \dontrun{}
```
