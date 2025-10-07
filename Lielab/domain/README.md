# Lielab::domain

Domain is the primary submodule and purpose of Lielab.

Class method styles:

- `get_`: Gets a value or property inherent to the mathematical structure or how Lielab handles the associated data.
- `set_`: Sets a value or property inherent to the mathematical structure or how Lielab handles the associated data.
- `from_`: Statically initializes the object from a convenient representation.
- `to_`: Returns values in a convenient representation.

## Lie Algebras

Lie algebras follow math conventions of being lower case. These maps are guaranteed:

- `x.get_vector()`: Returns the row vector representation. Real.
- `x.set_vector(v)`: Sets data from the row vector representation. Does not reshape `x`. Real.
- `x.get_matrix()`: Returns the matrix vector representation. Real or complex.
- One index call, `x(ind)`: Gets a value from the row vector representation. Always real.
- Two index call, `x(ind1, ind2)`: Gets a value from the square matrix representation. Real or complex.
- Math ops `+` and `-`: Operations on the vector space.

These maps are not guaranteed, but may exist depending on the object:

- `x.from_vector()`: Creates a new object of minimum shape from a row vector representation. Real.
- One index bracket, `x[ind]`: Returns a value given the context of the object.

## Lie Groups

Lie groups follow math conventions of being upper case. These maps are guaranteed:

- `g.serialize()`: Returns an unambiguous row data representation. Real.
- `g.unserialize(v)`: Sets data from the unambiguous row data representation. Does not reshape `g`. Real.
- `g.get_matrix()`: Returns the matrix representation. Real or complex.
- Two index call, `g(ind1, ind2)`: Gets a value from the square matrix representation. Real or complex.
- Math ops `*` and `g.inverse()`: Operations on the vector space.

These maps are not guaranteed, but may exist depending on the object:

- One index call, `g(ind)`: Returns a value given the context of the object.
- One index bracket, `x[ind]`: Returns a value given the context of the object.

## Smooth Manifolds

tbd...

## LieIII

Template macro providing compile time type information of the Lie algebra-group relationship, $\mathfrak{Lie}\mathrm{III} : \mathfrak{g} \leftrightarrow G$.

## Composite structures

CompositeX objects create composite structures preserving both mathematical properties and implementation details.

- `CompositeAlgebra`: Creates an object preserving algebra operations, $\mathfrak{g} = \mathfrak{h}_1 \oplus \mathfrak{h}_2 \oplus \cdots$.
- `CompositeGroup`: Creates an object preserving group operations, $G = H_1 \times H_2 \times \cdots$.
- `CompositeManifold`: Creates an object preserving smooth manifold operations, $M = N_1 \times N_2 \times \cdots$.
