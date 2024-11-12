# InvariantManifolds.jl

The purpose of this package is to provide a convenience tool to numerically investigate the low dimensional
invariant manifolds. We don't provide any numerical stability or reliable commitment. 

The main idea of this package is quite simple. By using the local manifolds of saddles, extend this manifolds
step by step. We just keep the points near enough (with the distance and curvature control) to ensure the accuracy of the numerical manifolds.

The two-dimensional algorithm for smooth mapping of this package is based on the one-dimensional algorithm. The numerical manifolds are represented as plenty near enough circles.

The most interesting part of this package is the computing of the non-smooth invariant manifolds in a reliable way. See the examples for more details.