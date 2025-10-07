#ifndef LIELAB_INTEGRATE_IVPMETHODS_RUNGEKUTTA_HPP
#define LIELAB_INTEGRATE_IVPMETHODS_RUNGEKUTTA_HPP

#include "../Coefficients.hpp"
#include "../ODESolution.hpp"

#include "IVPSettings.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/functions.hpp"

#include <Eigen/Core>

#include <functional>
#include <limits>
#include <memory>
#include <string>

namespace Lielab::integrate
{

enum class RungeKuttaStatus
{
    SUCCESS = 0,
    DO_STEP0 = 1,
    DO_STEP1 = 2,
    ESTIMATE_ERROR = 3,
};

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
    int status = 0;
    std::string message = "";

    double reltol = 1.0e-4;
    double abstol = 1.0e-8;
    double error_estimate;

    bool can_variable_step = false;
    size_t order = 0;


    // RK Params
    Eigen::MatrixXd A;
    Eigen::VectorXd B;
    Eigen::VectorXd Bhat;
    Eigen::VectorXd C;
    Eigen::VectorXd e;
    size_t n;
    
    // Other
    size_t stage = 0;
    double t0;
    double dt;
    std::vector<Eigen::VectorXd> K;

    double next_t;
    Eigen::VectorXd next_theta;
    Eigen::VectorXd next_theta2;

    RungeKutta();
    RungeKutta(const Coefficients tableau);

    virtual RungeKuttaStatus init(const double t0, const double dt, const Eigen::VectorXd& dy);
    virtual RungeKuttaStatus step_0();
    virtual RungeKuttaStatus step_1(const Eigen::VectorXd& dy);
    virtual RungeKuttaStatus postprocess();
    virtual RungeKuttaStatus estimate_error(const Eigen::VectorXd &yl, const Eigen::VectorXd &yh);
    Eigen::VectorXd operator()(const EuclideanIVP_vectorfield_t vf, const Eigen::VectorXd& y0, const double t0_, const double dt_);
};

enum class RungeKuttaFlowStatus
{
    ERROR_SMALL_DT = -5,
    ERROR_NEGATIVE_DT = -4,
    ERROR_SUCCEEDED_BUT_TOL_THO = -3,
    ERROR_MAX_ITERATIONS = -2,
    ERROR_INPUT = -1,
    SUCCESS = 0,
    DO_STEP0 = 1,
    DO_STEP1 = 2,
};

/*
* Numerically calculates the integral curve of a vector field through y0 from t0 to tf.
*/
class RungeKuttaFlow
{
    public:
    double small = 0.2;
    double large = 10.0;
    double pessimist = 0.9;
    double dt_min = 1e-4;
    double dt_max = 10.0;

    double reltol = 1e-4;
    double abstol = 1e-8;

    double dt = std::numeric_limits<double>::quiet_NaN();

    bool tolerance_not_met = false;
    bool variable_time_step = true;
    bool has_event = false;
    // bool rebase_every_step = true;

    double event_current = std::numeric_limits<double>::quiet_NaN();
    double event_next = std::numeric_limits<double>::quiet_NaN();

    std::shared_ptr<RungeKutta> method = std::make_shared<RungeKutta>(RungeKutta());
    Lielab::utils::newton search;

    int new_exact = 1;
    int num_step = 8;

    int iterations = 0;
    int max_iterations = -1;
    size_t tind;
    Eigen::VectorXd _tspan;
    double _error_estimate;

    double dt_recommend = std::numeric_limits<double>::quiet_NaN();
    double dt_save = std::numeric_limits<double>::quiet_NaN();

    double _tcurrent;
    Eigen::VectorXd _ycurrent;
    Eigen::VectorXd _ynext;

    Eigen::VectorXd _yl;
    Eigen::VectorXd _yh;
    double _err = std::numeric_limits<double>::quiet_NaN();

    // Output variables
    ODESolution _out;

    RungeKuttaFlowStatus init(const Eigen::VectorXd& tspan, const Eigen::VectorXd& y0);

    // bool step0(const Eigen::VectorXd& next_low);
    RungeKuttaFlowStatus step0(const Eigen::VectorXd& next_low, const double next_error);
    RungeKuttaFlowStatus step1();
    void stepE(double event_val);
    void postprocess();

    ODESolution operator()(const EuclideanIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Eigen::VectorXd& y0, const Lielab::integrate::IVPOptions& options);
};

}

#endif
