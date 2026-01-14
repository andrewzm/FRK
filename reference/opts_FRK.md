# FRK options

The main options list for the FRK package.

## Usage

``` r
opts_FRK
```

## Format

List of 2

- `$` :

  `set:function(opt,value)`

- `$` :

  `get:function(opt)`

## Details

`opts_FRK` is a list containing two functions, `set` and `get`, which
can be used to set options and retrieve options, respectively. Currently
`FRK` uses three options:

- "progress"::

  a flag indicating whether progress bars should be displayed or not

- "verbose"::

  a flag indicating whether certain progress messages should be shown or
  not. Currently this is the only option applicable to `method` = "TMB"

- "parallel"::

  an integer indicating the number of cores to use. A number 0 or 1
  indicates no parallelism

## Examples

``` r
opts_FRK$set("progress",1L)
opts_FRK$get("parallel")
#> [1] 0
```
