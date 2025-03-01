# InvariantManifolds.jl

This package aims to provide a convenient tool for numerically investigating low-dimensional invariant manifolds. We offer no guarantees regarding numerical stability or reliability.

**Main idea:**

The core concept of this package is straightforward. By utilizing the local manifolds of saddle points, we progressively extend these manifolds. We maintain points in close proximity (controlled by distance and curvature) to ensure the accuracy of the numerical manifolds.

In this package, the two-dimensional algorithm for smooth mappings is built upon the one-dimensional algorithm. Numerical manifolds are represented as a collection of sufficiently close points, forming circles.

The most compelling aspect of this package lies in its ability to reliably compute non-smooth invariant manifolds. Please refer to the examples for further details.