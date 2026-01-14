# Plotting themes

Formats a ggplot object for neat plotting.

## Usage

``` r
LinePlotTheme()

EmptyTheme()
```

## Value

Object of class `ggplot`

## Details

`LinePlotTheme()` creates `ggplot` object with a white background, a
relatively large font, and grid lines. `EmptyTheme()` on the other hand
creates a `ggplot` object with no axes or legends.

## Examples

``` r
if (FALSE) { # \dontrun{
X <- data.frame(x=runif(100),y = runif(100), z = runif(100))
LinePlotTheme() + geom_point(data=X,aes(x,y,colour=z))
EmptyTheme() + geom_point(data=X,aes(x,y,colour=z))} # }
```
