# InvariantManifolds.jl

[![](https://img.shields.io/badge/docs-online-blue.svg)](https://Xiaomingzzhang.github.io/InvariantManifolds.jl/dev/)
[![Build Status](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml?query=branch%3Amaster)

This package is designed to compute the invariant manifolds of dynamical systems.

**Main features:**

  - Computing one-dimensional manifolds of saddles for smooth mappings;
  - Computing two-dimensional manifolds of saddles for autonomous vector fields and smooth mappings;
  - Computing one-dimensional manifolds of saddles for non-smooth mappings: these mappings are time-T maps of non-smooth ODE systems, such as impact systems, piecewise smooth systems, and their combinations;
  - Computing two-dimensional manifolds of saddles for non-smooth autonomous vector fields.

To use this package, first install [Julia](https://julialang.org/), then run:

```julia
using Pkg;
Pkg.add(url="https://github.com/Xiaomingzzhang/InvariantManifolds.jl")
```

For more information, please refer to the [documentation](https://Xiaomingzzhang.github.io/InvariantManifolds.jl/dev/) of this package. The documentation was originally written in Chinese and then translated into English using AI texteditor Windsurf. A Chinese version of the documentation is also maintained for Chinese-speaking users. See [中文文档](https://xiaomingzzhang.github.io/cndocs/).