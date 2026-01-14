# MODIS cloud data

An image of a cloud taken by the Moderate Resolution Imaging
Spectroradiometer (MODIS) instrument aboard the Aqua satellite (MODIS
Characterization Support Team, 2015).

## Usage

``` r
MODIS_cloud_df
```

## Format

A data frame with 33,750 rows and 3 variables:

- x:

  x-coordinate

- y:

  y-coordinate

- z:

  binary dependent variable: 1 if cloud is present, 0 if no cloud. This
  variable has been thresholded from the original continuous measurement
  of radiance supplied by the MODIS instrument

- z_unthresholded:

  The original continuous measurement of radiance supplied by the MODIS
  instrument

## References

MODIS Characterization Support Team (2015). MODIS 500m Calibrated
Radiance Product.NASA MODIS Adaptive Processing System, Goddard Space
Flight Center, USA.
