#ifndef LIELAB_INTEGRATE_IVPMETHODS_SOLVE_IVP_HPP
#define LIELAB_INTEGRATE_IVPMETHODS_SOLVE_IVP_HPP

#include "IVPMethods/IVPCommon.hpp"

#include "Lielab/domain.hpp"

#include <Eigen/Core>

#include <functional>
#include <memory>
#include <limits>
#include <vector>

namespace Lielab::integrate
{

IVPSolution solve_ivp(const EuclideanIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Eigen::VectorXd& y0, const IVPOptions options = IVPOptions());
IVPSolution solve_ivp(const HomogeneousIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Lielab::domain::CompositeManifold& y0, const IVPOptions options = IVPOptions());

}

#endif
