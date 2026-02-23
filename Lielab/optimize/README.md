# Lielab::optimize

## Extremization

Main drivers:
- `minimize()`
- `maximize()`

Main methods:
- `LineSearch`: A very fast, but very innacurate line search method using only 1 jacobian function evaluation. Shouldn't be used on its own.
- `Newton`: A Newton method. Requires gradient.
- `Golden`: A 1-D Golden search. Kind of robust to NaNs.

## RootSolving

Main driver: `solve_roots()`

Main methods:
- `Newton`: A Newton method.
- `Hybrid`: A wrapper for Eigen's built in Hybrid solver.
