# World map

This world map was extracted from the package `maps` v.3.0.1 by running
`ggplot2::map_data("world")`. To reduce the data size, only every third
point of this data frame is contained in `worldmap`.

## Usage

``` r
worldmap
```

## Format

A data frame with 33971 rows and 6 variables:

- long:

  longitude coordinate

- lat:

  latitude coordinate

- group:

  polygon (region) number

- order:

  order of point in polygon boundary

- region:

  region name

- subregion:

  subregion name

## References

Original S code by Becker, R.A. and Wilks, R.A. This R version is by
Brownrigg, R. Enhancements have been made by Minka, T.P. and Deckmyn, A.
(2015) maps: Draw Geographical Maps, R package version 3.0.1.
