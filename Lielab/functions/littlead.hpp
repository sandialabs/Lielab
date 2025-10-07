#ifndef LIELAB_FUNCTIONS_LITTLEAD_HPP
#define LIELAB_FUNCTIONS_LITTLEAD_HPP

#include "commutator.hpp"

#include "Lielab/utils.hpp"

#include <vector>

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::glr ad_numerical(const LA & a, const int p = 1);

template <typename LA>
Lielab::domain::glr ad(const LA & a, const int p = 1);

template <typename LA>
LA ad_numerical(const LA & a, const LA & b, const int p = 1);

template <typename LA>
LA ad(const LA & a, const LA & b, const int p = 1);

template <>
Lielab::domain::glr ad(const Lielab::domain::cn & a, const int p);

template <>
Lielab::domain::cn ad(const Lielab::domain::cn & a, const Lielab::domain::cn & b, const int p);

template <>
Lielab::domain::glr ad(const Lielab::domain::glr & a, const int p);

template <>
Lielab::domain::glr ad(const Lielab::domain::glr & a, const Lielab::domain::glr & b, const int p);

// template <>
// Lielab::domain::glr ad(const Lielab::domain::glc & a, const int p);

// template <>
// Lielab::domain::glc ad(const Lielab::domain::glc & a, const Lielab::domain::glc & b, const int p);

template <>
Lielab::domain::glr ad(const Lielab::domain::rn & a, const int p);

template <>
Lielab::domain::rn ad(const Lielab::domain::rn & a, const Lielab::domain::rn & b, const int p);

template <>
Lielab::domain::glr ad(const Lielab::domain::so & a, const int p);

template <>
Lielab::domain::so ad(const Lielab::domain::so & a, const Lielab::domain::so & b, const int p);

template <>
Lielab::domain::glr ad(const Lielab::domain::su & a, const int p);

template <>
Lielab::domain::su ad(const Lielab::domain::su & a, const Lielab::domain::su & b, const int p);

template <>
Lielab::domain::glr ad(const Lielab::domain::se & a, const int p);

template <>
Lielab::domain::se ad(const Lielab::domain::se & a, const Lielab::domain::se & b, const int p);

}

#include "littlead.tpp"

#endif
