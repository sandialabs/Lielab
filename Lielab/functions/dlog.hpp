#ifndef LIELAB_FUNCTIONS_DLOG_HPP
#define LIELAB_FUNCTIONS_DLOG_HPP

#include "dexpinv.hpp"

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::glr dlog_numerical(const LA & a, const size_t order = 5);

template <typename LA>
Lielab::domain::glr dlog(const LA & a, const size_t order = 5);

template <typename LA>
LA dlog_numerical(const LA & a, const LA & b, const size_t order = 5);

template <typename LA>
LA dlog(const LA & a, const LA & b, const size_t order = 5);

}

#include "dlog.tpp"

#endif
