# Lielab::testing

Various checks and assertions for data.

## check_topology

A real-time check for topological equivalence between two manifolds.

``` python
>>> X = CompositeManifold([SU.identity(5), SO.identity(4)])
>>> Y = CompositeManifold([exp(su.basis(0, 5)), exp(so.basis(2, 4))])
>>> check_topology(X, Y)
True
>>> Z = CompositeManifold([SU(2), SO(3)])
>>> check_topology(X, Z)
False
```

## check_almost_equal_tol

A real-time check for equivalence of data to a specified tolerance.

``` python
>>> X = CompositeManifold([exp((1.0)*so.basis(0,3))])
>>> Y = CompositeManifold([exp((1.0 + 1e-14)*so.basis(0,3))])
>>> Z = CompositeManifold([exp((1.0 + 1e-13)*so.basis(0,3))])
>>> check_almost_equal_tol(X, Y)
True
>>> check_almost_equal_tol(X, Z)
False
>>> check_almost_equal_tol(X, Y, reltol=1e-15, abstol=1e-15)
False
```

Always runs check_topology first

``` python
>>> X = CompositeManifold([exp((1.0)*so.basis(0,3))])
>>> Y = CompositeManifold([exp((1.0 + 1e-14)*so.basis(0,3)), SU(2)])
>>> check_almost_equal_tol(X, Y)
False
```

## check_almost_equal_nulp

A real-time check for equivalence of data to a specified number of units of last place. Uses the larger of the
two numbers being checked as the scale for a single ULP.

``` python
>>> X = CompositeManifold([exp((1.0)*so.basis(0,3))])
>>> Y = CompositeManifold([exp((1.0 + 1e-14)*so.basis(0,3))])
>>> check_almost_equal_nulp(X, Y)
False
>>> check_almost_equal_nulp(X, Y, nulp=100)
True
```

Setting `gate=True` will use the larger of the two numbers *and* 1.0, so the smallest single ULP encountered
will typically be about `1e-16`. Good for checking quantities that could be near zero.

``` python
>>> X = CompositeAlgebra([so.from_vector([0, 0, 0])])
>>> Y = CompositeAlgebra([so.from_vector([1e-16, 0, 0])])
>>> check_almost_equal_nulp(X, Y)
False
>>> check_almost_equal_nulp(X, Y, gate=True)
True
```

Always runs check_topology first

``` python
>>> X = CompositeAlgebra([so.from_vector([0, 0, 0])])
>>> Y = CompositeAlgebra([so.from_vector([1e-16, 0, 0]), su(2)])
>>> check_almost_equal_nulp(X, Y, gate=True)
False
```

## lielab_assert

A preprocessor macro for cassert-like real-time checking of data that cannot be
verified at compile-time. Should be used to check if a problem was set
up incorrectly:

``` cpp
double myfunc(const Eigen::VectorXd& x, const Eigen::VectorXd& y)
{
    lielab_assert(x.size() == y.size(), "Size of x and y must be the same");
    ...
}
```

Should **not** be used in failure of a formula or algorithm to converge:

``` cpp
double mysolver::myfunc(const Eigen::VectorXd& x, const Eigen::VectorXd& y)
{
    lielab_assert(x.size() == y.size(), "Size of x and y must be the same");

    Eigen::VectorXd z(5);

    if (y(0) == 0.0)
    {
        // Handle divide by zero
        z(0) = std::numeric_limits<double>::quiet_NaN();
        this->message = "Algorithm failed: divide by zero.";
        return z; // Return some value and continue code execution.
    }
    ...
}
```

Unlike cassert, lielab_assert will provide context information about the failure

``` python
>>> solver.myfunc([1,2],[3,4,5])
RuntimeError    Traceback (most recent call last)
----> 1 solver.myfunc([1,2],[3,4,5])

RuntimeError: Lielab runtime error
Pass condition: x.size() == y.size()
Reason: Size of x and y must be the same
Where: myfunc() located in mysolver.cpp (Line 50)
```

lielab_assert is included when the `LIELAB_INCLUDE_ASSERTS` preprocessor variable is defined and will be present in all debug and release builds.
