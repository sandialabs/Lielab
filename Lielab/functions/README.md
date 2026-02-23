# Lielab::functions

Functions submodule of Lielab. Functions are in general guaranteed to work on all domains by the
numerical procedures. Some functions are overloaded for certain domains.

## exp()

Exponential function.

``` python
>>> Y = exp(so(3))
>>> Y.to_string()
SO(3)
>>> Z = exp(su(5))
>>> Z.to_string()
SU(5)
>>> G = exp(CompositeAlgebra([so(3), su(2), sp(4)]))
>>> G.to_string()
SO(3) x SU(2) x SP(4, R)
```

## log()

Logarithm function.

``` python
>>> y = log(SO(3))
>>> y.to_string()
so(3)
>>> z = log(SU(5))
>>> z.to_string()
su(5)
>>> x = log(CompositeGroup([SO(3), SU(2), SP(4)]))
>>> x.to_string()
so(3) ⊕ su(2) ⊕ sp(4, R)
```
