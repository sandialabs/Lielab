#ifndef LIELAB_FUNCTIONS_COMMUTATOR_HPP
#define LIELAB_FUNCTIONS_COMMUTATOR_HPP

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <typename LA>
LA commutator(const LA & a, const LA & b);

template <>
Lielab::domain::CompositeAlgebra commutator(const Lielab::domain::CompositeAlgebra &a, const Lielab::domain::CompositeAlgebra &b);

}

#include "commutator.tpp"

#endif
