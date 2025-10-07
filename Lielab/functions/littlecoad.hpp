#ifndef LIELAB_FUNCTIONS_LITTLECOAD_HPP
#define LIELAB_FUNCTIONS_LITTLECOAD_HPP

#include "commutator.hpp"
#include "littlead.hpp"

#include "Lielab/utils.hpp"

#include <vector>

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::glr coad_numerical(const LA & x, const int p = 1);

template <typename LA>
LA coad_numerical(const LA & x, const LA & y, const int p = 1);

template <typename LA>
Lielab::domain::glr coad(const LA & x, const int p = 1);

template <typename LA>
LA coad(const LA & x, const LA & y, const int p = 1);

}

#include "littlecoad.tpp"

#endif
