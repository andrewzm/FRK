# Draw a map of the world with country boundaries.

Layers a `ggplot2` map of the world over the current `ggplot2` object.

## Usage

``` r
draw_world(g = ggplot() + theme_bw() + xlab("") + ylab(""), inc_border = TRUE)
```

## Arguments

- g:

  initial ggplot object

- inc_border:

  flag indicating whether a map border should be drawn or not; see
  details.

## Details

This function uses
[`ggplot2::map_data()`](https://ggplot2.tidyverse.org/reference/map_data.html)
in order to create a world map. Since, by default, this creates lines
crossing the world at the (-180,180) longitude boundary, the function
`.homogenise_maps()` is used to split the polygons at this boundary into
two. If `inc_border` is TRUE, then a border is drawn around the lon-lat
space; this option is most useful for projections that do not yield
rectangular plots (e.g., the sinusoidal global projection).

## See also

the help file for the dataset [`worldmap`](worldmap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
draw_world(g = ggplot())} # }
```
