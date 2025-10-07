#include "RungeKutta.hpp"

#include "../ODESolution.hpp"
#include "IVPSettings.hpp"

#include <Eigen/Core>

#include <tuple>

namespace Lielab::integrate
{

// IVPMethod::IVPMethod()
// {
//     /*!
//     * Creates a Method.
//     */
    
// }

// Eigen::VectorXd IVPMethod::operator()(const std::function<Eigen::VectorXd(const double, const Eigen::VectorXd&)> vectorfield, const Eigen::VectorXd& y0, const double t0, const double dt)
// {
//     /*!
//      * Primary usage of IVPMethod. Do not call this.
//      */

//     this->error_estimate = std::numeric_limits<double>::quiet_NaN();
//     return y0;
// }

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
* Author / Date: Sparapany / 2020
*/


// RungeKutta::RungeKutta() : Method()
// {
//     /*! \f{equation*}{ () \rightarrow RungeKutta \f}
//     * Instantiates a new RungeKutta object.
//     */

// }

RungeKutta::RungeKutta()
{
    /*! \f{equation*}{ () \rightarrow MuntheKaas \f}
    * Instantiates a new MuntheKaas object.
    */
    
    const auto [_A, _b, _bhat, _c, _e, _order, _stages, _variable] = get_butcher_tableau(Coefficients::RKV87r);
    
    this->A = _A;
    this->B = _b;
    this->Bhat = _bhat;
    this->C = _c;
    this->e = _e;
    this->n = _stages;
    this->order = _order;
    this->can_variable_step = _variable;
}

RungeKutta::RungeKutta(const Coefficients tableau)
{
    /*!
    * Instantiates a new RungeKutta object with a specific RK method.
    */
    
    const auto [_A, _b, _bhat, _c, _e, _order, _stages, _variable] = get_butcher_tableau(tableau);
    
    this->A = _A;
    this->B = _b;
    this->Bhat = _bhat;
    this->C = _c;
    this->e = _e;
    this->n = _stages;
    this->order = _order;
    this->can_variable_step = _variable;
}

RungeKuttaStatus RungeKutta::init(const double t0_, const double dt_, const Eigen::VectorXd& dy)
{
    /*! \f{equation*}{ (\mathbb{R}, M, \mathbb{R}, \mathfrak{g}^{n}) \rightarrow () \f}
    * Initializes the RungeKaas process.
    */

    using namespace Lielab::domain;

    this->stage = 0;

    this->t0 = t0_;
    this->dt = dt_;

    this->K = std::vector<Eigen::VectorXd>(this->n);

    this->next_t = t0;
    this->next_theta = Eigen::VectorXd::Zero(dy.size());

    return RungeKuttaStatus::DO_STEP0;
}

RungeKuttaStatus RungeKutta::step_0()
{
    /*!
    * Advances the solution to current stage and returns pair
    * (next_t, next_y) to be evaluated by the vectorfield.
    */
    
    if (this->stage == 0)
    {
        this->next_t = this->t0;
        return RungeKuttaStatus::DO_STEP1;
    }

    this->next_t = this->t0 + this->dt*this->C(this->stage);

    this->next_theta *= 0.0;
    for (size_t jj = 0; jj < this->stage; jj++)
    {
        this->next_theta += this->dt*this->A(this->stage, jj)*this->K[jj];
    }

    return RungeKuttaStatus::DO_STEP1;
}

RungeKuttaStatus RungeKutta::step_1(const Eigen::VectorXd& dy)
{
    /*!
    * Advance solution
    */

    this->K[this->stage] = dy;

    if (this->stage == 0)
    {
        // Initialize U if this was the first stage.
        this->next_theta = Eigen::VectorXd::Zero(this->K[0].size());
    }

    this->stage++;

    if (this->stage == n)
    {
        // Finished integration
        return RungeKuttaStatus::SUCCESS;
    }

    // Needs another step
    return RungeKuttaStatus::DO_STEP0;
}

RungeKuttaStatus RungeKutta::postprocess()
{
    using namespace Lielab::domain;

    this->next_theta *= 0.0;

    for (size_t ii = 0; ii < this->n; ii++)
    {
        this->next_theta += this->dt*this->B(ii)*this->K[ii];
    }

    this->error_estimate = 0.0;

    if (this->can_variable_step == true)
    {
        const size_t sz = this->next_theta.size();
        this->next_theta2 = Eigen::VectorXd::Zero(sz);
        for (size_t ii = 0; ii < this->n; ii++)
        {
            this->next_theta2 += this->dt*this->Bhat(ii)*this->K[ii];
        }
        return RungeKuttaStatus::ESTIMATE_ERROR;
    }

    return RungeKuttaStatus::SUCCESS;
}

RungeKuttaStatus RungeKutta::estimate_error(const Eigen::VectorXd &y1, const Eigen::VectorXd &y2)
{
    const Eigen::VectorXd scale = this->abstol + this->reltol*y1.array().max(y2.array()).cwiseAbs();
    const size_t sz = y1.size();

    Eigen::VectorXd err0 = Eigen::VectorXd::Zero(sz);
    for (size_t ii = 0; ii < n; ii++)
    {
        err0 += this->K[ii]*this->e[ii]*this->dt;
    }

    const Eigen::VectorXd err = err0.array()/scale.array();
    this->error_estimate = err.norm()/std::sqrt(static_cast<double>(sz));

    return RungeKuttaStatus::SUCCESS;
}

Eigen::VectorXd RungeKutta::operator()(const EuclideanIVP_vectorfield_t vf, const Eigen::VectorXd& y0, const double t0_, const double dt_)
{
    /*!
     * Primary usage of RungeKutta.
     */

    this->status = 0;
    this->message = "";

    Eigen::VectorXd dy = vf(t0_, y0);
    Eigen::VectorXd next_y = 0.0*y0;
    RungeKuttaStatus status = this->init(t0_, dt_, dy);
    
    while (static_cast<int>(status) > 0)
    {
        if (status == RungeKuttaStatus::DO_STEP0)
        {
            status = this->step_0();
        }
        else if (status == RungeKuttaStatus::DO_STEP1)
        {
            next_y = y0 + this->next_theta;
            dy = vf(this->next_t, next_y);
            if (dy.array().isNaN().any())
            {
                this->status = -1;
                this->message = "NaNs in vectorfield.";
                return y0;
            }
            else if (dy.array().isInf().any())
            {
                this->status = -2;
                this->message = "Infs in vectorfield.";
                return y0;
            }

            status = this->step_1(dy);
        }
    }

    status = this->postprocess();
    next_y = y0 + this->next_theta;

    if (status == RungeKuttaStatus::ESTIMATE_ERROR)
    {
        status = this->estimate_error(y0 + this->next_theta, y0 + this->next_theta2);
    }

    return next_y;
}

/*
* The Flow class. Numerically calculates the integral curve of a vector field through y0 from t0 to tf.
*/

RungeKuttaFlowStatus RungeKuttaFlow::init(const Eigen::VectorXd& tspan, const Eigen::VectorXd& y0)
{
    if (tspan.size() < 2)
    {
        return RungeKuttaFlowStatus::ERROR_INPUT;
    }

    if ((this->variable_time_step == true) && (this->method->can_variable_step == false))
    {
        this->variable_time_step = false;
    }

    if (std::isnan(this->dt))
    {
        this->dt = 0.01*(tspan(1) - tspan(0));
    }

    this->dt = std::max(this->dt_min, std::min(this->dt_max, this->dt));
    this->_tspan = tspan;
    this->tind = 1;
    const size_t num_eoms = y0.size();
    
    this->_tcurrent = tspan(0);
    this->_ycurrent = y0;

    this->iterations = 0;

    this->_out = ODESolution(num_eoms);
    this->_out.add_data(tspan(0), y0);

    if (this->dt <= 0.0)
    {
        return RungeKuttaFlowStatus::ERROR_NEGATIVE_DT;
    }

    return RungeKuttaFlowStatus::DO_STEP0;
}

RungeKuttaFlowStatus RungeKuttaFlow::step0(const Eigen::VectorXd& next_low, const double next_error)
{
    /*
     * Checks if we will accept the current approximation by method. If not, comes up with a new one.
     */

    const double eps = std::numeric_limits<double>::epsilon();

    if (this->variable_time_step)
    {
        double dt_new_mult = 1.0;

        if (next_error < 1.0)
        {
            // Step accepted

            if (std::abs(next_error) < 4*eps)
            {
                // Prevent 0 to a power if error is 0.
                dt_new_mult = this->large;
            }
            else
            {
                dt_new_mult = std::min(this->large, this->pessimist*std::pow(next_error, -1.0/(static_cast<double>(method->order) + 1.0)));
            }

            this->dt_recommend = dt_new_mult*this->dt;
            
            // Obey min and max dt amounts
            this->dt_recommend = std::min(this->dt_recommend, this->dt_max);
            this->dt_recommend = std::max(this->dt_recommend, this->dt_min);
            this->_ynext = next_low;
            return RungeKuttaFlowStatus::DO_STEP1;
        }
        else
        {
            // Step rejected

            if (std::abs(this->dt - this->dt_min) < 4.0*eps)
            {
                // Accept step: dt is at lower bound. Continue solution as normal but warn.
                this->_ynext = next_low;
                this->tolerance_not_met = true;
                return RungeKuttaFlowStatus::DO_STEP1;
            }

            dt_new_mult = std::max(this->small, this->pessimist*std::pow(next_error, -1.0/(static_cast<double>(method->order) + 1.0)));
            this->dt = dt_new_mult*this->dt;
            
            // Obey min and max dt amounts
            this->dt = std::min(this->dt, this->dt_max);
            this->dt = std::max(this->dt, this->dt_min);
            return RungeKuttaFlowStatus::DO_STEP0;
        }
    }

    this->_ynext = next_low;
    if (next_error >= 1.0)
    {
        // Accept step but signal warning (despite being fixed step).
        this->tolerance_not_met = false;
        return RungeKuttaFlowStatus::DO_STEP1;
    }

    // Accept step.
    return RungeKuttaFlowStatus::DO_STEP1;
}

RungeKuttaFlowStatus RungeKuttaFlow::step1()
{
    /*
     * Solution was accepted. Save.
     */

    const double eps = std::numeric_limits<double>::epsilon();

    this->iterations += 1;
    this->_tcurrent += this->dt;
    Eigen::VectorXd ydiff = this->_ynext;
    this->_ycurrent = this->_ynext;
    this->_out.add_data(this->_tcurrent, this->_ycurrent, ydiff);

    if (!std::isnan(this->dt_save))
    {
        // If a saved timestep is stored, reset it.
        this->dt = this->dt_save;
        this->dt_save = std::numeric_limits<double>::quiet_NaN();
    }
    else if (!std::isnan(this->dt_recommend))
    {
        // If a new dt was recommended, use it.
        this->dt = this->dt_recommend;
        this->dt_recommend = std::numeric_limits<double>::quiet_NaN();
    }

    // Advance tind if there are more in tspan
    if (this->_tcurrent > this->_tspan(tind) - 1e-10)
    {
        // If we're passed the last value in tspan, we're done.
        if (this->_tcurrent > (this->_tspan.tail<1>()(0) - 1e-10))
        {
            return RungeKuttaFlowStatus::SUCCESS;
        }

        tind += 1;
    }

    // Error check. dt should be within bounds by this point.
    // TODO: This should never get thrown, yet here we are.
    if (this->dt < 100*eps)
    {
        return RungeKuttaFlowStatus::ERROR_SMALL_DT;
    }

    if (!std::isnan(this->event_next) && this->event_next <= 0)
    {
        return RungeKuttaFlowStatus::SUCCESS;
    }

    // Check if the next time step will cross the next time in tspan
    if (this->_tcurrent + this->dt > (this->_tspan(tind) - 1e-10))
    {
        // Save the current dt
        this->dt_save = this->dt;

        // Set the next dt to end exactly on the next tspan value
        this->dt = this->_tspan(tind) - this->_tcurrent;
    }

    if (this->max_iterations > 0 && iterations >= this->max_iterations)
    {
        // Max iterations exceeded.
        return RungeKuttaFlowStatus::ERROR_MAX_ITERATIONS;
    }

    // Nothing was caught, keep running.
    return RungeKuttaFlowStatus::DO_STEP0;
}

void RungeKuttaFlow::stepE(double event_val)
{
    // if (event_val < 0 || std::abs(event_val) < tol)
    // {
    //     int status = 0; // TODO: Return this
    // }
}

void RungeKuttaFlow::postprocess()
{
    _out.trim_chunk(iterations);
}

ODESolution RungeKuttaFlow::operator()(const EuclideanIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Eigen::VectorXd& y0, const IVPOptions & options)
{
    /*!
    * The main evaluation method for the Flow class.
    */

    // Check if the solution we're running has an event
    const double event_val = dynamics.event(tspan(0), y0);
    this->has_event = false;
    if (!std::isnan(event_val))
    {
        this->has_event = true;
    }

    // Initialize the Flow object
    RungeKuttaFlowStatus status = this->init(tspan, y0);

    // Main loop. Run until algorithm says it's done
    while (static_cast<int>(status) > 0)
    {
        if (status == RungeKuttaFlowStatus::DO_STEP0)
        {
            Eigen::VectorXd ynext = this->method->operator()(dynamics.vectorfield, this->_ycurrent, this->_tcurrent, this->dt);
            
            if (this->has_event)
            {
                this->event_current = dynamics.event(this->_tcurrent, this->_ycurrent);
                this->event_next = dynamics.event(this->_tcurrent + this->dt, ynext);
                if ((this->event_current >= 0) && (this->event_next <= 0))
                {
                    // Event crossed
                    this->search.lower = std::numeric_limits<double>::epsilon();
                    this->search.upper = this->dt;
                    const auto fun = [&](const double _dt)
                    {
                        const Eigen::VectorXd _ynext = this->method->operator()(dynamics.vectorfield, this->_ycurrent, this->_tcurrent, _dt);
                        return dynamics.event(this->_tcurrent + _dt, _ynext);
                    };

                    const double xopt = this->search(fun, this->dt/2.0);
                    this->dt = xopt;
                    ynext = this->method->operator()(dynamics.vectorfield, this->_ycurrent, this->_tcurrent, this->dt);
                    this->event_next = dynamics.event(this->_tcurrent + this->dt, ynext);
                }
            }

            status = this->step0(ynext, this->method->error_estimate);
        }
        else if (status == RungeKuttaFlowStatus::DO_STEP1)
        {
            status = this->step1();
        }
    }

    this->postprocess();

    if (status == RungeKuttaFlowStatus::SUCCESS && this->tolerance_not_met)
    {
        status = RungeKuttaFlowStatus::ERROR_SUCCEEDED_BUT_TOL_THO;
    }

    ODESolution out = ODESolution(this->_out);
    out.status = static_cast<int>(status);
    if (status == RungeKuttaFlowStatus::SUCCESS)
    {
        out.message = "Integration succeeded.";
    }
    else if (status == RungeKuttaFlowStatus::ERROR_INPUT)
    {
        out.message = "Input error: tspan must have size greater than 2.";
    }
    else if (status == RungeKuttaFlowStatus::ERROR_MAX_ITERATIONS)
    {
        out.message = "Max iterations exceeded.";
    }
    else if (status == RungeKuttaFlowStatus::ERROR_SUCCEEDED_BUT_TOL_THO)
    {
        out.message = "Integration succeeded, but tolerance not met.";
    }
    else if (status == RungeKuttaFlowStatus::ERROR_NEGATIVE_DT)
    {
        out.message = "dt 0 or negative. Only positive dt is supported.";
    }
    else if (status == RungeKuttaFlowStatus::ERROR_SMALL_DT)
    {
        out.message = "dt became too small.";
    }

    return out;
}

}
