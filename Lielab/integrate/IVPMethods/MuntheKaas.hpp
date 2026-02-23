#ifndef LIELAB_INTEGRATE_IVPMETHODS_MUNTHEKAAS_HPP
#define LIELAB_INTEGRATE_IVPMETHODS_MUNTHEKAAS_HPP

#include "Coefficients.hpp"
#include "IVPCommon.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/functions.hpp"
#include "Lielab/optimize.hpp"

#include <Eigen/Core>

#include <functional>
#include <limits>
#include <memory>
#include <string>

namespace Lielab::integrate
{

class MuntheKaas
{
    public:
    IVPStatus status = IVPStatus::RUNNING;
    bool success;
    std::string message = "";

    double abstol = 1.0e-8;
    double reltol = 1.0e-4;
    double error_estimate = std::numeric_limits<double>::quiet_NaN();

    bool can_variable_step = false;
    bool implicit = false;
    int order = 0;


    // RK Params
    Eigen::MatrixXd A;
    Eigen::VectorXd B;
    Eigen::VectorXd Bhat;
    Eigen::VectorXd C;
    Eigen::VectorXd e;
    int n;
    
    // Other
    Eigen::MatrixXd K;

    MuntheKaas(const RungeKuttaCoefficients tableau = RungeKuttaCoefficients::RKV87r);
    void estimate_error(const Eigen::VectorXd& yl, const Eigen::VectorXd& yh, const double dt);
    Lielab::domain::CompositeManifold operator()(const HomogeneousIVPSystem& dynamics, const Lielab::domain::CompositeManifold& y0, const double t0, const double dt);
};

/*
* Numerically calculates the integral curve of a vector field through y0 from t0 to tf.
*/
class MuntheKaasFlow
{
    public:
    std::string message = "";
    IVPStatus status = IVPStatus::RUNNING;
    bool success = false;
    
    bool tolerance_not_met = false;
    int iterations = 0;
    
    Eigen::VectorXd _ynext;

    // Output variables
    IVPSolution operator()(const HomogeneousIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Lielab::domain::CompositeManifold& y0, const IVPOptions& options);
};

}

#endif
