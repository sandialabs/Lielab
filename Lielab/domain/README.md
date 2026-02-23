# Lielab::domain

Domain is the primary submodule and purpose of Lielab. The following objects and methods are implemented to a mixed degree of testing and reliability.

## Smooth Manifolds

Every object is a smooth manifold. Implemented operations include:

### Empty constructor

Creates an empty default object with

``` python
>>> y = SE()
<lielab.domain.SE>
```

### .to_string()

Returns a mathematical representation as a string.

``` python
>>> y = SE()
>>> print(y.to_string())
SE(0)
```

### .get_dimension()

Returns the dimension of the manifold.

``` python
>>> y = SE(2)
>>> y.get_dimension()
3
>>> z = SU(5)
>>> z.get_dimension()
24
```

### .get_size()

Returns the size of the array returned by `.serialize()` and expected by `.unserialize(...)`.

``` python
>>> y = SE(2)
>>> y.get_size()
6
>>> y.serialize()
[0., 0., 1., 0., 0., 1.]
```

### .get_point()

Returns the point on the manifold in its internal representation.

``` python
>>> y = SE(3)
>>> y.to_string()
SE(3)
>>> point = manifold.get_point()
(<lielab.domain.RN>, <lielab.domain.SO>)
>>> print(point[0].to_string() + ' x ' + point[1].to_string())
R^3 x SO(3)
```

### .serialize()

Returns an unambiguous representation of the data point as an array of real numbers.

``` python
>>> y = SO(3)
>>> y.serialize()
[1., 0., 0., 0., 1., 0., 0., 0., 1.]
```

This operation ignores manifold constraints.

### .unserialize()

Directly assigns data to the manifold point using the unambiguous real array representation.

``` python
>>> y = SO(3)
>>> y.unserialize([1,2,3,4,5,6,7,8,9])
>>> y.serialize()
[1., 2., 3., 4., 5., 6., 7., 8., 9.]
```

This operation ignores manifold constraints.

## Lie Algebras

Lie algebras are generally lowercase. Implemented operations include

### .is_abelian()

Returns whether or not a Lie algebra element is Abelian.

``` python
>>> y = rn(3)
>>> y.is_abelian()
True
>>> z = se(3)
>>> z.is_abelian()
False
>>> a = CompositeAlgebra([rn(3)])
>>> a.is_abelian()
True
>>> b = CompositeAlgebra([rn(3), se(3)])
>>> b.is_abelian()
False
```

### Matrix constructor

Allows a Lie algebra to be constructed from a matrix.

``` python
>>> y = se(np.random.uniform(size=(4,4)))
>>> y.get_dimension()
6
>>> y.get_matrix()
[[0.03509347, 0.86758351, 0.49301856, 0.50478953],
 [0.55413825, 0.66343378, 0.15511134, 0.29276696],
 [0.29083431, 0.046435  , 0.87869363, 0.37540926],
 [0.        , 0.        , 0.        , 0.        ]]
```

The input matrix will not be projected and may violate manifold constraints.

### .basis()

Instantiate a Lie algebra from a specific basis element with signature

$$\textrm{basis}(\textrm{index}, \textrm{shape})$$

``` python
>>> y = so.basis(0,3)
>>> y.get_matrix()
[[ 0., -0.,  0.],
 [ 0.,  0., -1.],
 [-0.,  1.,  0.]]
```

### .zero()

Instantiate a Lie algebra at the zero element by shape

$$\textrm{zero}(\textrm{shape})$$

``` python
>>> y = so.zero(4)
>>> y.get_matrix()
[[0., 0., 0., 0.],
 [0., 0., 0., 0.],
 [0., 0., 0., 0.],
 [0., 0., 0., 0.]]
```

### .from_vector()

Instantiate a Lie algebra from a vector. Will internally choose a proper shape to fit the data.

``` python
>>> y = su.from_vector([1,2,3])
>>> y.get_shape()
2
>>> y.get_matrix()
[[ 0.+3.j, -2.+1.j],
 [ 2.+1.j,  0.-3.j]]
```

### .get_shape()

Returns the size of the Lie algebra in its matrix representation such that $\textrm{get\_matrix} : \mathfrak{g} \rightarrow \mathbb{F}^{\textrm{shape} \; \times \; \textrm{shape}}$

``` python
>>> y = su.from_vector([1,2,3,4,5,6,7,8])
>>> y.get_shape()
3
>>> y.get_matrix()
[[ 0.+11.61880215j, -4. +1.j        , -5. +2.j        ],
 [ 4. +1.j        ,  0. -2.38119785j, -6. +3.j        ],
 [ 5. +2.j        ,  6. +3.j        ,  0. -9.23760431j]]
```

### .get_matrix()

Returns the matrix (Ado's) representation in a structure that preserves the commutator using matrix operations.

``` python
>>> x = so.from_vector([1, 0, 0])
>>> y = so.from_vector([0, 1, 0])
>>> zhat = np.dot(x.get_matrix(), y.get_matrix()) - np.dot(y.get_matrix(), x.get_matrix())
[[ 0., -1.,  0.],
 [ 1.,  0.,  0.],
 [ 0.,  0.,  0.]]
>>> so(zhat).get_vector()
[-0., 0., 1.]
```

### .get_vector()

Returns a mininum dimensional real row vector representation of the data.

``` python
>>> x = so.basis(0,4)
>>> x.get_dimension()
6
>>> x.get_vector()
[1., 0., 0., 0., 0., 0.]
```

### .set_vector()

Assigns data using the minimum dimensional real row vector representation.

``` python
>>> y = su(2)
>>> y.get_dimension()
3
>>> y.set_vector([0, 0, 1])
>>> y.get_matrix()
[[0.+1.j, 0.+0.j],
 [0.+0.j, 0.-1.j]]
```

### Lie algebra math operations

Addition, subtraction, and scalar multiplication and division are all implemented.

``` python
>>> x = so.basis(0, 3)
>>> y = so.basis(1, 3)
>>> z = 2.5*x + y/3.0 + (-x)
>>> z.get_vector()
[1.5, 0.33333333, -0.]
```

## Lie Groups

Lie groups are generally uppercase. Implemented operations include:

### .is_abelian()

Returns whether or not a Lie algebra element is Abelian.

``` python
>>> y = RN(3)
>>> y.is_abelian()
True
>>> z = SE(3)
>>> z.is_abelian()
False
>>> a = CompositeGroup([RN(3)])
>>> a.is_abelian()
True
>>> b = CompositeGroup([RN(3), SE(3)])
>>> b.is_abelian()
False
```

### Matrix constructor

Allows a Lie algebra to be constructed from a matrix.

``` python
>>> y = SE(np.random.uniform(size=(4,4)))
>>> y.get_dimension()
6
>>> y.get_matrix()
[[0.66769104, 0.49394359, 0.52606881, 0.59960518],
 [0.47506925, 0.4294394 , 0.47229971, 0.22036006],
 [0.2788966 , 0.8844491 , 0.51738232, 0.76277922],
 [0.        , 0.        , 0.        , 1.        ]]
```

The input matrix will not be projected and may violate manifold constraints.

### .identity()

Instantiate a Lie group at the identity element by shape

$$\textrm{identity}(\textrm{shape})$$

``` python
>>> y = SO.identity(4)
>>> y.get_matrix()
[[1., 0., 0., 0.],
 [0., 1., 0., 0.],
 [0., 0., 1., 0.],
 [0., 0., 0., 1.]]
```

### .get_shape()

Returns the size of the Lie group in its matrix representation such that $\textrm{get\_matrix} : G \rightarrow \mathbb{F}^{\textrm{shape} \; \times \; \textrm{shape}}$

``` python
>>> y = SU(3))
>>> y.get_shape()
3
>>> y.get_matrix()
[[1.+0.j, 0.+0.j, 0.+0.j],
 [0.+0.j, 1.+0.j, 0.+0.j],
 [0.+0.j, 0.+0.j, 1.+0.j]]
```

### .get_matrix()

Returns the matrix (Ado's) representation in a structure that preserves the group using matrix operations.

``` python
>>> X = exp(so.from_vector([1, 0, 0]))
>>> X.to_string()
SO(3)
>>> X.inverse().get_matrix()
[[ 1.        ,  0.        ,  0.        ],
 [ 0.        ,  0.54030231,  0.84147098],
 [ 0.        , -0.84147098,  0.54030231]]
>>> np.linalg.inv(X.get_matrix())
[[ 1.        ,  0.        ,  0.        ],
 [ 0.        ,  0.54030231,  0.84147098],
 [-0.        , -0.84147098,  0.54030231]]
```

### Lie group operations

Product and inverse operations are implemented for Lie groups.

``` python
>>> Y = exp(se.from_vector([1.0, 0.0, 0.0, 2.0, 0.0, 0.0]))
>>> (Y*Y*Y).get_matrix()
[[ 1.        ,  0.        ,  0.        ,  3.        ],
 [ 0.        ,  0.96017029,  0.2794155 ,  0.        ],
 [ 0.        , -0.2794155 ,  0.96017029,  0.        ],
 [ 0.        ,  0.        ,  0.        ,  1.        ]]
>>> (Y*Y.inverse()).get_matrix()
[[1., 0., 0., 0.],
 [0., 1., 0., 0.],
 [0., 0., 1., 0.],
 [0., 0., 0., 1.]]
```

## Composite objects

Composite objects are direct products of elements giving an array-like interface. Implemented objects include:

### CompositeAlgebra

A container preserving the aforementioned Lie algebra operations, $\mathfrak{g} = \mathfrak{h}_1 \oplus \mathfrak{h}_2 \oplus \cdots$

### CompositeGroup

A container preserving the aforementioned Lie group operations, $G = H_1 \times H_2 \times \cdots$

### CompositeManifold

A container preserving the aforementioned Smooth Manifold operations, $M = N_1 \times N_2 \times \cdots$

## Additional methods

Some objects will have additional undocumented helper methods given their context. For example with SO(3):

``` python
>>> R = exp(so.from_vector([np.pi/8, np.pi/8, 0]))
>>> R.to_eulerangles_body123()
[0.9595312915239355, 0.6861731241754524, -0.36753179727574403]
>>> R.to_eulerangles_space123()
[0.41347024105782615, 0.3820481415498507, 0.08107219749699343]
>>> R.to_eulerangles_body121()
[0.19889203328807747, 0.3901415209103397, 0.19889203328807747]
>>> R.to_quaternion()
[0.9616939461100744, 0.1938359538569061, 0.1938359538569061, 0.0]
```

I try to add errors where possible

``` python
>>> SO(5).to_quaternion()
RuntimeError    Traceback (most recent call last)
----> 1 SO(5).to_quaternion()

RuntimeError: Lielab runtime error
Pass condition: this->get_shape() == 3
Reason: Expected shape 3. Got 5.
Where: to_quaternion() located in SO.cpp (Line X)
```

## LieIII (C++ Only)

Template macro providing compile time type information of the Lie algebra-group relationship, $\mathfrak{Lie}\mathrm{III} : \mathfrak{g} \leftrightarrow G$.
