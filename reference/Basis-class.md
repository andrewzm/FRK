# Basis functions

An object of class `Basis` contains the basis functions used to
construct the matrix \\S\\ in FRK.

## Details

Basis functions are a central component of `FRK`, and the package is
designed to work with user-defined specifications of these. For
convenience, however, several functions are available to aid the user to
construct a basis set for a given set of data points. Please see
[`auto_basis`](auto_basis.md) for more details. The function
[`local_basis`](local_basis.md) helps the user construct a set of local
basis functions (e.g., bisquare functions) from a collection of location
and scale parameters.

## Slots

- `manifold`:

  an object of class `manifold` that contains information on the
  manifold and the distance measure used on the manifold. See
  [`manifold-class`](manifold-class.md) for more details

- `n`:

  the number of basis functions in this set

- `fn`:

  a list of length `n`, with each item the function of a specific basis
  function

- `pars`:

  a list of parameters where the \\i\\-th item in the list contains the
  parameters of the \\i\\-th basis function, `fn[[i]]`

- `df`:

  a data frame containing other attributes specific to each basis
  function (for example the geometric centre of the local basis
  function)

- `regular`:

  logical indicating if the basis functions (of each resolution) are in
  a regular grid

## See also

[`auto_basis`](auto_basis.md) for automatically constructing basis
functions and [`show_basis`](show_basis.md) for visualising basis
functions.
