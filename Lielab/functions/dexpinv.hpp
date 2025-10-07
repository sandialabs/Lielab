#ifndef LIELAB_FUNCTIONS_DEXPINV_HPP
#define LIELAB_FUNCTIONS_DEXPINV_HPP

#include "dexp.hpp"
#include "littlead.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::glr dexpinv_numerical(const LA & a, const size_t order = 5);

template <typename LA>
Lielab::domain::glr dexpinv(const LA & a, const size_t order = 5);

template <typename LA>
LA dexpinv_numerical(const LA & a, const LA & b, const size_t order = 5);

template <typename LA>
LA dexpinv(const LA & a, const LA & b, const size_t order = 5);

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::glr & x, const size_t order);

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::glr & x, const Lielab::domain::glr & y, const size_t order);

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::so & x, const size_t order);

template <>
Lielab::domain::so dexpinv(const Lielab::domain::so & a, const Lielab::domain::so & b, const size_t order);

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::se & y, const size_t order);

template <>
Lielab::domain::se dexpinv(const Lielab::domain::se & a, const Lielab::domain::se & b, const size_t order);

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::su & x, const size_t order);

template <>
Lielab::domain::su dexpinv(const Lielab::domain::su & a, const Lielab::domain::su & b, const size_t order);

template <>
Lielab::domain::CompositeAlgebra dexpinv(const Lielab::domain::CompositeAlgebra & a, const Lielab::domain::CompositeAlgebra & b, const size_t order);

}

#include "dexpinv.tpp"

#endif
