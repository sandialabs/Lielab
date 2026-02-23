#ifndef LIELAB_INTEGRATE_IVPMETHODS_RUNGEKUTTA_HPP
#define LIELAB_INTEGRATE_IVPMETHODS_RUNGEKUTTA_HPP

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

/*!
* A Runge-Kutta type Method.
*
* Let \f$M\f$ be a homogeneous manifold with action \f$L : (G, M) \rightarrow M\f$,
* let \f$t \in \mathbb{R}\f$ be time, and let \f$f : (t, M) \rightarrow \mathfrak{g}\f$ be a set of first-order ordinary
* differential equations. A Munthe-Kaas method solves the
* differential equations by solving the related set of equations:
* 
* \f{equation*}{
* \frac{dU}{dt} = d\phi_U^{-1}(f(t, y))
* \f}
* 
* where \f$ \phi : \mathfrak{g} \rightarrow G \f$ is a coordinate map and
* \f$ d\phi_U^{-1} : \mathfrak{g} \rightarrow T \mathfrak{g} \f$. Given an \f$s\f$-stage
* Runge-Kutta method, the iterative Munthe-Kaas method is summarized as
* 
* \f{eqnarray*}{
* U_i &=& dt \sum_{j=1}^{i-1} a_{ij} K_j \\
* y_i &=& L(\phi(U_i), y_{i-1}) \\
* K_i &=& d\phi_U^{-1} f(t, y_i) \\
* i &=& 1, \dots, s \\
* V &=& dt \sum_{i=1}^s b_i K_i \\
* y_f &=& L(\phi(V), y_0)
* \f}
* 
* Author / Date: Sparapany / May 2025
*/

class RungeKutta
{
    public:
    std::string message = "";
    IVPStatus status = IVPStatus::RUNNING;
    bool success;

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

    RungeKutta(const RungeKuttaCoefficients tableau = RungeKuttaCoefficients::RKV87r);
    void estimate_error(const Eigen::VectorXd& yl, const Eigen::VectorXd& yh, const double dt);
    Eigen::VectorXd operator()(const EuclideanIVP_vectorfield_t vf, const Eigen::VectorXd& y0, const double t0, const double dt);
};

/*
* Numerically calculates the integral curve of a vector field through y0 from t0 to tf.
*/
class RungeKuttaFlow
{
    public:
    std::string message = "";
    IVPStatus status = IVPStatus::RUNNING;
    bool success = false;

    bool tolerance_not_met = false;
    int iterations = 0;
    
    Eigen::VectorXd _ynext;

    // Output variables
    IVPSolution operator()(const EuclideanIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Eigen::VectorXd& y0, const IVPOptions& options);
};

}

#endif
