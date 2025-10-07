#ifndef LIELAB_FUNCTIONS_CAYLEY_HPP
#define LIELAB_FUNCTIONS_CAYLEY_HPP

#include "Lielab/domain.hpp"

namespace Lielab::functions
{
template <typename LA>
Lielab::domain::LieIII<LA> Cayley(const LA & g);

template <>
Lielab::domain::CN Cayley(const Lielab::domain::cn & a);

template <>
Lielab::domain::GLC Cayley(const Lielab::domain::glc & a);

template <>
Lielab::domain::SU Cayley(const Lielab::domain::su & a);

Lielab::domain::CompositeGroup Cayley(const Lielab::domain::CompositeAlgebra & la);

template <typename LA>
Lielab::domain::LieIII<LA> Cayley2(const LA & g);

template <typename LA>
LA dCayleyinv(const LA & u, const LA & v);

Lielab::domain::CompositeAlgebra dCayleyinv(const Lielab::domain::CompositeAlgebra & a, const Lielab::domain::CompositeAlgebra & b);

}

#include "Cayley.tpp"

#endif
