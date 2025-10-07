#ifndef LIELAB_FUNCTIONS_EXP_HPP
#define LIELAB_FUNCTIONS_EXP_HPP

#include "Lielab/domain.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::functions
{

template<typename LA>
Lielab::domain::LieIII<LA> exp_numerical(const LA & la);

template<typename LA>
Lielab::domain::LieIII<LA> exp(const LA & la);

template<>
Lielab::domain::CN exp(const Lielab::domain::cn & la);

template<>
Lielab::domain::GLR exp(const Lielab::domain::glr & la);

template<>
Lielab::domain::GLC exp(const Lielab::domain::glc & la);

template<>
Lielab::domain::RN exp(const Lielab::domain::rn & la);

template<>
Lielab::domain::SO exp(const Lielab::domain::so & x);

template<>
Lielab::domain::SE exp(const Lielab::domain::se & y);

template<>
Lielab::domain::CompositeGroup exp(const Lielab::domain::CompositeAlgebra & la);

}

#include "exp.tpp"

#endif
