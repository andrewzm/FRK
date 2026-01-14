# manifold

The class `manifold` is virtual; other manifold classes inherit from
this class.

## Details

A `manifold` object is characterised by a character variable `type`,
which contains a description of the manifold, and a variable `measure`
of type `measure`. A typical measure is the Euclidean distance.

`FRK` supports five manifolds; the real line (in one dimension),
instantiated by using [`real_line()`](real_line.md); the 2D plane,
instantiated by using [`plane()`](plane.md); the 2D-sphere surface S2,
instantiated by using [`sphere()`](sphere.md); the R2 space-time
manifold, instantiated by using [`STplane()`](STplane.md), and the S2
space-time manifold, instantiated by using [`STsphere()`](STsphere.md).
User-specific manifolds can also be specified, however helper functions
that are manifold specific, such as `auto_BAUs` and `auto_basis`, only
work with the pre-configured manifolds. Importantly, one can change the
distance function used on the manifold to synthesise anisotropy or
heterogeneity. See the vignette for one such example.

## See also

[`real_line`](real_line.md), [`plane`](plane.md), [`sphere`](sphere.md),
[`STplane`](STplane.md) and [`STsphere`](STsphere.md) for constructing
manifolds.
