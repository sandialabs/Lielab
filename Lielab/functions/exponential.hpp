#ifndef LIELAB_FUNCTIONS_EXPONENTIAL_HPP
#define LIELAB_FUNCTIONS_EXPONENTIAL_HPP

#include "Lielab/domain.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::functions
{

// exp
template <typename g> Lielab::domain::LieIII<g> exp_numerical(const g& x);
template <typename g> Lielab::domain::LieIII<g> exp(const g& x);

template <> Lielab::domain::CN exp(const Lielab::domain::cn& x);
template <> Lielab::domain::GLC exp(const Lielab::domain::glc& x);
template <> Lielab::domain::GLR exp(const Lielab::domain::glr& x);
template <> Lielab::domain::RN exp(const Lielab::domain::rn& x);
template <> Lielab::domain::SE exp(const Lielab::domain::se& x);
template <> Lielab::domain::SO exp(const Lielab::domain::so& x);
// sp
// su

// log
template <typename G> Lielab::domain::LieIII<G> log_numerical(const G& g);
template <typename G> Lielab::domain::LieIII<G> log(const G& g);

template <> Lielab::domain::cn log(const Lielab::domain::CN& g);
template <> Lielab::domain::rn log(const Lielab::domain::RN& g);
template <> Lielab::domain::so log(const Lielab::domain::SO& g);
template <> Lielab::domain::se log(const Lielab::domain::SE& g);

// dexp
template <typename g, typename out_t = Lielab::domain::glr> out_t dexp_numerical(const g& x, const int order = 5);
template <typename g, typename out_t = Lielab::domain::glr> out_t dexp(const g& x, const int order = 5);
template <typename g> g dexp_numerical(const g& x, const g& y, const int order = 5);
template <typename g> g dexp(const g& x, const g& y, const int order = 5);

template <> Lielab::domain::glr dexp(const Lielab::domain::glr& x, const int order);
template <> Lielab::domain::glr dexp(const Lielab::domain::glr& x, const Lielab::domain::glr& y, const int order);
template <> Lielab::domain::glr dexp(const Lielab::domain::se& x, const int order);
template <> Lielab::domain::se dexp(const Lielab::domain::se& x, const Lielab::domain::se& y, const int order);
template <> Lielab::domain::glr dexp(const Lielab::domain::so& x, const int order);
template <> Lielab::domain::so dexp(const Lielab::domain::so& x, const Lielab::domain::so& y, const int order);
template <> Lielab::domain::glr dexp(const Lielab::domain::su& x, const int order);
template <> Lielab::domain::su dexp(const Lielab::domain::su& x, const Lielab::domain::su& y, const int order);

// dexpinv
template <typename g, typename out_t = Lielab::domain::glr> out_t dexpinv_numerical(const g& x, const int order = 5);
template <typename g, typename out_t = Lielab::domain::glr> out_t dexpinv(const g& x, const int order = 5);
template <typename g> g dexpinv_numerical(const g& x, const g& y, const int order = 5);
template <typename g> g dexpinv(const g& x, const g& y, const int order = 5);

template <> Lielab::domain::glr dexpinv(const Lielab::domain::glr& x, const int order);
template <> Lielab::domain::glr dexpinv(const Lielab::domain::glr& x, const Lielab::domain::glr& y, const int order);
template <> Lielab::domain::glr dexpinv(const Lielab::domain::se& x, const int order);
template <> Lielab::domain::se dexpinv(const Lielab::domain::se& x, const Lielab::domain::se& y, const int order);
template <> Lielab::domain::glr dexpinv(const Lielab::domain::so& x, const int order);
template <> Lielab::domain::so dexpinv(const Lielab::domain::so& x, const Lielab::domain::so& y, const int order);
template <> Lielab::domain::glr dexpinv(const Lielab::domain::su& x, const int order);
template <> Lielab::domain::su dexpinv(const Lielab::domain::su& x, const Lielab::domain::su& y, const int order);

// dlog
template <typename g, typename out_t = Lielab::domain::glr> out_t dlog_numerical(const g& x, const int order = 5);
template <typename g, typename out_t = Lielab::domain::glr> out_t dlog(const g& x, const int order = 5);
template <typename g> g dlog_numerical(const g& x, const g& y, const int order = 5);
template <typename g> g dlog(const g& x, const g& y, const int order = 5);


// dloginv
template <typename g, typename out_t = Lielab::domain::glr> out_t dloginv_numerical(const g& x, const int order = 5);
template <typename g, typename out_t = Lielab::domain::glr> out_t dloginv(const g& x, const int order = 5);
template <typename g> g dloginv_numerical(const g& x, const g& y, const int order = 5);
template <typename g> g dloginv(const g& x, const g& y, const int order = 5);

// Composite overloads
template <> Lielab::domain::CompositeGroup exp_numerical(const Lielab::domain::CompositeAlgebra& x);
template <> Lielab::domain::CompositeGroup exp(const Lielab::domain::CompositeAlgebra& x);
template <> Lielab::domain::CompositeAlgebra log_numerical(const Lielab::domain::CompositeGroup& g);
template <> Lielab::domain::CompositeAlgebra log(const Lielab::domain::CompositeGroup& g);
template <> Lielab::domain::CompositeAlgebra dexp_numerical(const Lielab::domain::CompositeAlgebra& x, const int order);
template <> Lielab::domain::CompositeAlgebra dexp(const Lielab::domain::CompositeAlgebra& x, const int order);
template <> Lielab::domain::CompositeAlgebra dexp_numerical(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order);
template <> Lielab::domain::CompositeAlgebra dexp(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order);
template <> Lielab::domain::CompositeAlgebra dexpinv_numerical(const Lielab::domain::CompositeAlgebra& x, const int order);
template <> Lielab::domain::CompositeAlgebra dexpinv(const Lielab::domain::CompositeAlgebra& x, const int order);
template <> Lielab::domain::CompositeAlgebra dexpinv_numerical(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order);
template <> Lielab::domain::CompositeAlgebra dexpinv(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order);
template <> Lielab::domain::CompositeAlgebra dlog_numerical(const Lielab::domain::CompositeAlgebra& x, const int order);
template <> Lielab::domain::CompositeAlgebra dlog(const Lielab::domain::CompositeAlgebra& x, const int order);
template <> Lielab::domain::CompositeAlgebra dlog_numerical(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order);
template <> Lielab::domain::CompositeAlgebra dlog(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order);
template <> Lielab::domain::CompositeAlgebra dloginv_numerical(const Lielab::domain::CompositeAlgebra& x, const int order);
template <> Lielab::domain::CompositeAlgebra dloginv(const Lielab::domain::CompositeAlgebra& x, const int order);
template <> Lielab::domain::CompositeAlgebra dloginv_numerical(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order);
template <> Lielab::domain::CompositeAlgebra dloginv(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order);

}

#include "exponential.tpp"

#endif
