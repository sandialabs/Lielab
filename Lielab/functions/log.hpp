#ifndef LIELAB_FUNCTIONS_LOG_HPP
#define LIELAB_FUNCTIONS_LOG_HPP

#include "Lielab/domain.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::functions
{

template<typename LG>
Lielab::domain::LieIII<LG> log_numerical(const LG& G);

template<typename LG>
Lielab::domain::LieIII<LG> log(const LG& G);

template <>
Lielab::domain::cn log(const Lielab::domain::CN& G);

template <>
Lielab::domain::rn log(const Lielab::domain::RN& G);

template <>
Lielab::domain::so log(const Lielab::domain::SO& W);

template<>
Lielab::domain::se log(const Lielab::domain::SE& Y);

Lielab::domain::CompositeAlgebra log(const Lielab::domain::CompositeGroup& M);

}

#include "log.tpp"

#endif
