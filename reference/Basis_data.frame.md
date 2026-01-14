# Basis-function data frame object

Tools for retrieving and manipulating the data frame within Basis
objects. Use the assignment `data.frame()<-` with care; no checks are
made to ensure the data frame conforms with the object.

## Usage

``` r
data.frame(x) <- value

# S4 method for class 'Basis'
x$name

# S4 method for class 'Basis'
x$name <- value

# S4 method for class 'Basis'
data.frame(x) <- value

# S4 method for class 'TensorP_Basis'
data.frame(x) <- value

# S3 method for class 'Basis'
as.data.frame(x, ...)

# S3 method for class 'TensorP_Basis'
as.data.frame(x, ...)
```

## Arguments

- x:

  the obect of class `Basis` we are assigning the new data to or
  retrieving data from

- value:

  the new data being assigned to the Basis object

- name:

  the field name to which values will be retrieved or assigned inside
  the Basis object's data frame

- ...:

  unused

## Examples

``` r
G <- local_basis()
df <- data.frame(G)
print(df$res)
#> [1] 1
df$res <- 2
data.frame(G) <- df
```
