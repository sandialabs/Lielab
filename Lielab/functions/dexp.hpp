#ifndef LIELAB_FUNCTIONS_DEXP_HPP
#define LIELAB_FUNCTIONS_DEXP_HPP

#include "littlead.hpp"

#include "commutator.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/utils.hpp"

#include <cmath>

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::glr dexp_numerical(const LA & a, const size_t order = 5);

template <typename LA>
Lielab::domain::glr dexp(const LA & a, const size_t order = 5);

template <typename LA>
LA dexp_numerical(const LA & a, const LA & b, const size_t order = 5);

template <typename LA>
LA dexp(const LA & a, const LA & b, const size_t order = 5);

template <>
Lielab::domain::glr dexp(const Lielab::domain::glr & x, const size_t order);

template <>
Lielab::domain::glr dexp(const Lielab::domain::glr & x, const Lielab::domain::glr & y, const size_t order);

template <>
Lielab::domain::glr dexp(const Lielab::domain::so & x, const size_t order);

template <>
Lielab::domain::so dexp(const Lielab::domain::so & a, const Lielab::domain::so & b, const size_t order);

template <>
Lielab::domain::glr dexp(const Lielab::domain::se & y, const size_t order);

template <>
Lielab::domain::se dexp(const Lielab::domain::se & a, const Lielab::domain::se & b, const size_t order);

template <>
Lielab::domain::glr dexp(const Lielab::domain::su & x, const size_t order);

template <>
Lielab::domain::su dexp(const Lielab::domain::su & a, const Lielab::domain::su & b, const size_t order);

template <>
Lielab::domain::CompositeAlgebra dexp(const Lielab::domain::CompositeAlgebra & a, const Lielab::domain::CompositeAlgebra & b, const size_t order);

}

#include "dexp.tpp"

#endif
