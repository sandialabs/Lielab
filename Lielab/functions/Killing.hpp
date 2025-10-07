#ifndef LIELAB_FUNCTIONS_KILLING_HPP
#define LIELAB_FUNCTIONS_KILLING_HPP

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <typename LA>
double Killing(const LA & a, const LA & b);

template <typename LA>
Eigen::MatrixXd Killingform(const LA & g);

}

#include "Killing.tpp"

#endif
