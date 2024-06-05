# InvariantManifolds.jl

The purpose of this package is to provide a convenience tool to numerically investigate the low dimensional
invariant manifolds. We don't provide any numerical stability or reliable commitment. 

The main idea of this package is quite simple. By using the local manifolds of saddles, extend this manifolds
step by step. We just keep the points near enough to ensure the accuracy of the numerical manifolds. Usually,
the curvatures of nearby points are also considered in most algorithms. We only use the distance control since the curvatures requirement often can't be satisfied after several iterations, even in most simple example.
