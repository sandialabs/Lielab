#ifndef LIELAB_FUNCTIONS_AD_HPP
#define LIELAB_FUNCTIONS_AD_HPP

#include "Lielab/domain.hpp"

#include <cassert>

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::GLR Ad_numerical(const LA & a);

template <typename LA>
Lielab::domain::GLR Ad(const LA & a);

// template <typename LG>
// Lielab::domain::GLR Ad_numerical(const LG & A);

// template <typename LG>
// Lielab::domain::GLR Ad(const LG & A);

template <typename LA>
LA Ad(const LA & a, const LA & b);

template <typename LA>
LA Ad(const Lielab::domain::LieIII<LA> & A, const LA & b);

Lielab::domain::CompositeAlgebra Ad(const Lielab::domain::CompositeGroup & A, const Lielab::domain::CompositeAlgebra & b);

}

#include "Ad.tpp"

#endif
