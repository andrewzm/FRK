# ISEA Aperture 3 Hexagon (ISEA3H) Discrete Global Grid

The data used here were obtained from
https://webpages.sou.edu/~sahrk/dgg/isea.old/gen/isea3h.html and
represent ISEA discrete global grids (DGGRIDs) generated using the
`DGGRID` software. The original .gen files were converted to a data
frame using the function `dggrid_gen_to_df`, available with the
`dggrids` package. Only resolutions 0–6 are supplied with `FRK` and note
that resolution 0 of ISEA3H is equal to resolution 1 in `FRK`. For
higher resolutions `dggrids` can be installed from
`https://github.com/andrewzm/dggrids/` using `devtools`.

## Usage

``` r
isea3h
```

## Format

A data frame with 284,208 rows and 5 variables:

- id:

  grid identification number within the given resolution

- lon:

  longitude coordinate

- lat:

  latitude coordinate

- res:

  DGGRID resolution (0 – 6)

- centroid:

  A 0-1 variable, indicating whether the point describes the centroid of
  the polygon, or whether it is a boundary point of the polygon

## References

Sahr, K. (2008). Location coding on icosahedral aperture 3 hexagon
discrete global grids. Computers, Environment and Urban Systems, 32,
174–187.
