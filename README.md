# InvariantManifolds.jl

[![](https://img.shields.io/badge/docs-online-blue.svg)](https://Xiaomingzzhang.github.io/InvariantManifolds.jl/dev/)
[![Build Status](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml?query=branch%3Amaster)

This is a small package to compute the invariant manifolds of a dynamical system.

Main features:

- Compute saddles' one-dimensional manifolds of smooth mapping;
- Compute saddles' two-dimensional manifolds of autonomous vector field;
- Compute saddles' one-dimensional manifolds of non-smooth mapping: these mapping are the time-T-map of a non-smooth ODE systems, e.g., impact systems, piecewise smooth systems, and simple Fillippov systems.

To use this package, first install [julia](https://julialang.org/), then run 
```julia
using Pkg;
Pkg.add("https://github.com/Xiaomingzzhang/InvariantManifolds.jl")
```

For more information, see the [docs](https://Xiaomingzzhang.github.io/InvariantManifolds.jl/dev/) of this package.