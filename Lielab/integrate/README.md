# Lielab::integrate

Integration submodule of Lielab.

## IVPMethods

Methods for solving initial value problems.

Types of IVP problems:

- `EuclideanIVPSystem`: An IVP evolving on the classical Euclidean space $\mathbb{R}^n$. Roughly translates to the usual definion style with NumPy and Eigen.
- `HomogeneousIVPSystem`: An IVP evolving on a homogeneous manifold $M$. Requires Lielab::domain objects.
