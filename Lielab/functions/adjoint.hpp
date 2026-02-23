#ifndef LIELAB_FUNCTIONS_ADJOINT_HPP
#define LIELAB_FUNCTIONS_ADJOINT_HPP

#include "Lielab/domain.hpp"

#include <vector>

namespace Lielab::functions
{

// commutator
template <typename g> g commutator(const g& x, const g& y);

// ad
template <typename g, typename out_t = Lielab::domain::glr> out_t ad_numerical(const g& x, const int power = 1);
template <typename g, typename out_t = Lielab::domain::glr> out_t ad(const g& x, const int power = 1);
template <typename g> g ad_numerical(const g& x, const g& y, const int power = 1);
template <typename g> g ad(const g& x, const g& y, const int power = 1);

template <> Lielab::domain::glr ad(const Lielab::domain::cn& x, const int power);
template <> Lielab::domain::cn ad(const Lielab::domain::cn& x, const Lielab::domain::cn& y, const int power);
// template <> Lielab::domain::glr ad(const Lielab::domain::glc& x, const int power);
// template <> Lielab::domain::glc ad(const Lielab::domain::glc& x, const Lielab::domain::glc& y, const int power);
template <> Lielab::domain::glr ad(const Lielab::domain::glr& x, const int power);
template <> Lielab::domain::glr ad(const Lielab::domain::glr& x, const Lielab::domain::glr& y, const int power);
template <> Lielab::domain::glr ad(const Lielab::domain::rn& x, const int power);
template <> Lielab::domain::rn ad(const Lielab::domain::rn& x, const Lielab::domain::rn& y, const int power);
template <> Lielab::domain::glr ad(const Lielab::domain::se& x, const int power);
template <> Lielab::domain::se ad(const Lielab::domain::se& x, const Lielab::domain::se& y, const int power);
template <> Lielab::domain::glr ad(const Lielab::domain::so& x, const int power);
template <> Lielab::domain::so ad(const Lielab::domain::so& x, const Lielab::domain::so& y, const int power);
template <> Lielab::domain::glr ad(const Lielab::domain::su& x, const int power);
template <> Lielab::domain::su ad(const Lielab::domain::su& x, const Lielab::domain::su& y, const int power);

// Ad
// template <typename LA> Lielab::domain::GLR Ad_numerical(const LA& x);
// template <typename LA> Lielab::domain::GLR Ad(const LA& x);
// template <typename LG> Lielab::domain::GLR Ad_numerical(const LG& g);
// template <typename LG> Lielab::domain::GLR Ad(const LG& g);

// template <typename LA> LA Ad(const LA& a, const LA& b);

template <typename G> Lielab::domain::LieIII<G> Ad(const G& g, const Lielab::domain::LieIII<G>& x);

// coad
template <typename g, typename out_t = Lielab::domain::glr> out_t coad_numerical(const g& x, const int power = 1);
template <typename g, typename out_t = Lielab::domain::glr> out_t coad(const g& x, const int power = 1);
template <typename g> g coad_numerical(const g& x, const g& y, const int power = 1);
template <typename g> g coad(const g& x, const g& y, const int power = 1);

// Composite overloads
template <> Lielab::domain::CompositeAlgebra commutator(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y);

template <> Lielab::domain::CompositeAlgebra ad_numerical(const Lielab::domain::CompositeAlgebra& x, const int power);
template <> Lielab::domain::CompositeAlgebra ad_numerical(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int power);
template <> Lielab::domain::CompositeAlgebra ad(const Lielab::domain::CompositeAlgebra& x, const int power);
template <> Lielab::domain::CompositeAlgebra ad(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int power);

template <> Lielab::domain::CompositeAlgebra coad_numerical(const Lielab::domain::CompositeAlgebra& x, const int power);
template <> Lielab::domain::CompositeAlgebra coad_numerical(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int power);
template <> Lielab::domain::CompositeAlgebra coad(const Lielab::domain::CompositeAlgebra& x, const int power);
template <> Lielab::domain::CompositeAlgebra coad(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int power);

template <> Lielab::domain::CompositeAlgebra Ad(const Lielab::domain::CompositeGroup& g, const Lielab::domain::CompositeAlgebra& x);
}

#include "adjoint.tpp"

#endif
