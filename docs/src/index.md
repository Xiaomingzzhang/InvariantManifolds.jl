# InvariantManifolds.jl

This package provides tools for numerically investigating low-dimensional invariant
manifolds.

**Main idea:**

The core concept of this package is straightforward. By utilizing the local manifolds of saddle points, we progressively extend these manifolds. We maintain points in close proximity (controlled by distance and curvature) to ensure the accuracy of the numerical manifolds.

In this package, the two-dimensional algorithm for smooth mappings is built upon the one-dimensional algorithm. Numerical manifolds are represented as a collection of sufficiently close points, forming circles.

The continuation uses a scale-aware Goodman--Wróbel angle-times-length control,
followed by postprocessing thinning by default.  For non-smooth systems the
derivative of the time map is computed with saltation (jump) matrices, so saddle
periodic orbits that collide with impact surfaces can be located accurately.  See
the example pages and [Curve thinning and saltation matrices](extensions.md) for
details.
