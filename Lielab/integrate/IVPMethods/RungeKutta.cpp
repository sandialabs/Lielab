#include "RungeKutta.hpp"

#include "IVPCommon.hpp"

#include "Lielab/optimize.hpp"

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

RungeKutta::RungeKutta(const RungeKuttaCoefficients tableau)
{
    /*!
    * Instantiates a new RungeKutta object with a specific RK method.
    */
    
    const auto [_A, _b, _bhat, _c, _e, _order, _stages, _variable, _implicit] = get_butcher_tableau(tableau);
    
    this->A.resize(_stages, _stages);
    this->B.resize(_stages);
    this->Bhat.resize(_stages);
    this->C.resize(_stages);
    this->e.resize(_stages);
    this->n = _stages;
    this->order = _order;
    this->can_variable_step = _variable;
    this->implicit = _implicit;

    for (int ii = 0; ii < this->n; ii++)
    {
        for (int jj = 0; jj < this->n; jj++)
        {
            this->A(ii,jj) = _A[ii][jj];
        }
        this->B(ii) = _b[ii];
        this->Bhat(ii) = _bhat[ii];
        this->C(ii) = _c[ii];
        this->e(ii) = _e[ii];
    }
}

void RungeKutta::estimate_error(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2, const double dt)
{
    const Eigen::VectorXd scale = this->abstol + this->reltol*(y1.array().abs()).max(y2.array().abs());
    const int sz = static_cast<int>(y1.size());

    Eigen::VectorXd err0 = Eigen::VectorXd::Zero(sz);
    for (int ii = 0; ii < this->n; ii++)
    {
        err0 += this->K.row(ii)*this->e[ii]*dt;
    }

    const Eigen::VectorXd err = err0.array()/scale.array();
    this->error_estimate = err.norm()/std::sqrt(static_cast<double>(sz));
}

Eigen::VectorXd RungeKutta::operator()(const EuclideanIVP_vectorfield_t vf, const Eigen::VectorXd& y0, const double t0, const double dt)
{
    /*!
     * Primary usage of RungeKutta.
     */

    const Eigen::VectorXd dy0 = vf(t0, y0);
    this->success = false;

    if (dy0.array().isNaN().any())
    {
        this->status = IVPStatus::ERROR_NANS_IN_VF;
        this->message = "RungeKutta: NaNs in vectorfield.";
        this->success = false;
        return y0;
    }

    if (dy0.array().isInf().any())
    {
        this->status = IVPStatus::ERROR_INFS_IN_VF;
        this->message = "RungeKutta: Infs in vectorfield.";
        this->success = false;
        return y0;
    }

    const size_t n_eoms = y0.size();

    if (n_eoms == 0) return Eigen::VectorXd::Zero(0);

    this->K.noalias() = Eigen::MatrixXd::Zero(this->n, n_eoms);

    double next_t = t0;
    Eigen::VectorXd next_theta = Eigen::VectorXd::Zero(y0.size());

    int max_iter = 50;

    if (!this->implicit)
    {
        for (int stage = 0; stage < this->n; stage++)
        {
            if (stage == 0)
            {
                next_t = t0;
                next_theta.noalias() = y0*0.0;
            }
            else
            {
                next_t = t0 + dt*this->C(stage);
                next_theta *= 0.0;
                for (int jj = 0; jj < stage; jj++)
                {
                    next_theta += dt*this->A(stage, jj)*this->K.row(jj);
                }
            }

            this->K.row(stage).noalias() = vf(next_t, y0 + next_theta);
        }

        if (this->K.array().isNaN().any())
        {
            this->status = IVPStatus::ERROR_NANS_IN_VF;
            this->message = "RungeKutta: NaNs in vectorfield.";
            this->success = false;
            return y0;
        }

        if (this->K.array().isInf().any())
        {
            this->status = IVPStatus::ERROR_INFS_IN_VF;
            this->message = "RungeKutta: Infs in vectorfield.";
            this->success = false;
            return y0;
        }
    }
    else
    {
        int iteration = 0;
        double err = 9e9;
        const double tol = 1e-15;

        // for (size_t stage = 0; stage < this->n; stage++)
        // {
        //     this->K.row(stage) = dy0;
        // }

        while (err > tol && iteration < max_iter)
        {
            Eigen::MatrixXd Knew = this->K*1.0;
            for (int stage = 0; stage < this->n; stage++)
            {
                next_t = t0 + dt*this->C(stage);
                next_theta *= 0.0;
                for (int jj = 0; jj < this->n; jj++)
                {
                    next_theta += dt*this->A(stage, jj)*this->K.row(jj);
                }

                Knew.row(stage).noalias() = vf(next_t, y0 + next_theta);
            }

            err = std::sqrt((Knew - this->K).array().square().sum());
            this->K = Knew;

            if (this->K.array().isNaN().any())
            {
                this->status = IVPStatus::ERROR_NANS_IN_VF;
                this->message = "RungeKutta: NaNs in vectorfield.";
                this->success = false;
                return y0;
            }

            if (this->K.array().isInf().any())
            {
                this->status = IVPStatus::ERROR_INFS_IN_VF;
                this->message = "RungeKutta: Infs in vectorfield.";
                this->success = false;
                return y0;
            }

            iteration++;
        }
        
        // TODO: Error if max iter here?
    }


    next_theta *= 0.0;
    Eigen::VectorXd next_theta2 = 0.0*next_theta;
    for (int ii = 0; ii < this->n; ii++)
    {
        next_theta += dt*this->B(ii)*this->K.row(ii);
    }

    // Postprocess
    this->error_estimate = 0.0;
    if (this->can_variable_step == true)
    {
        const int sz = static_cast<int>(next_theta.size());
        next_theta2 = Eigen::VectorXd::Zero(sz);
        for (int ii = 0; ii < this->n; ii++)
        {
            next_theta2 += dt*this->Bhat(ii)*this->K.row(ii);
        }
        this->estimate_error(y0 + next_theta, y0 + next_theta2, dt);
    }

    this->status = IVPStatus::SUCCESS;
    this->success = true;
    this->message = "RungeKutta: Terminated successfully.";

    return y0 + next_theta;
}

/*
* The Flow class. Numerically calculates the integral curve of a vector field through y0 from t0 to tf.
*/

IVPSolution RungeKuttaFlow::operator()(const EuclideanIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Eigen::VectorXd& y0, const IVPOptions& options)
{
    /*!
    * The main evaluation method for the Flow class.
    */
    
    using namespace Lielab::optimize;

    // Preprocessing
    const ptrdiff_t n_eoms = y0.size();

    IVPSolution out(n_eoms);
    this->iterations = 0;
    this->message = "";
    this->success = false;
    this->status = IVPStatus::RUNNING;

    RungeKutta solver = RungeKutta(options.coefficients);

    const double dt_min = std::abs(options.dt_min);
    const double dt_max = std::abs(options.dt_max);

    solver.abstol = options.abstol;
    solver.reltol = options.reltol;
    double dt_save = std::numeric_limits<double>::quiet_NaN();
    double dt_recommend = std::numeric_limits<double>::quiet_NaN();

    // Check if the solution we're running has an event
    const double event_val = dynamics.event(tspan(0), y0);
    bool has_event = false;
    if (!std::isnan(event_val))
    {
        has_event = true;
    }

    double event_current = std::numeric_limits<double>::quiet_NaN();
    double event_next = std::numeric_limits<double>::quiet_NaN();

    // Initialize the Flow object
    int computestatus = 1;

    const double tf = tspan.tail<1>()(0);
    ptrdiff_t tind = 1;

    bool variable_time_step = options.variable_time_step;
    if (!solver.can_variable_step)
    {
        variable_time_step = false;
    }

    double dt = options.dt;
    if (std::isnan(dt))
    {
        dt = 0.01*(tspan(1) - tspan(0));
    }

    dt = std::max(dt_min, std::min(dt_max, dt));
    
    double t = tspan(0);
    Eigen::VectorXd y = y0;
    Eigen::VectorXd ynext = y;

    out.add_data(tspan(0), y);

    // Main loop. Run until algorithm says it's done
    while (computestatus == 1 || computestatus == 2)
    {
        if (computestatus == 1)
        {
            // Construct next possible step
            ynext.noalias() = solver.operator()(dynamics.vectorfield, y, t, dt);
            
            if (!solver.success)
            {
                // Error in method. Exit immediately and report the error.
                this->status = solver.status;
                this->message = solver.message;
                this->success = solver.success;
                out.trim_chunk();
                return out;
            }
            
            if (has_event)
            {
                event_current = dynamics.event(t, y);
                event_next = dynamics.event(t + dt, ynext);
                if ((event_current >= 0) && (event_next <= 0))
                {
                    // Event crossed. Root solve for 0 and exit immediately.
                    GoldenMinimize search;
                    ExtremizationOptions optimoptions;
                    optimoptions.abstol = options.abstol;
                    optimoptions.reltol = options.reltol;
                    optimoptions.max_iterations = 50;

                    const auto fun = [&](const Eigen::VectorXd& x)
                    {
                        const Eigen::VectorXd _ynext = solver.operator()(dynamics.vectorfield, y, t, x(0));
                        return std::abs(dynamics.event(t + x(0), _ynext));
                    };

                    EuclideanExtremizationSystem system(fun);
                    system.lower_bound = Lielab::utils::to_VectorXd({std::numeric_limits<double>::epsilon()});
                    system.upper_bound = Lielab::utils::to_VectorXd({dt});
                    const Eigen::VectorXd xopt = search(system, optimoptions);
                    dt = xopt(0);
                    
                    ynext.noalias() = solver.operator()(dynamics.vectorfield, y, t, dt);
                    
                    this->iterations += 1;
                    t += dt;
                    y.noalias() = ynext;
                    out.add_data(t, y);

                    // Exit
                    if (!this->tolerance_not_met)
                    {
                        this->status = IVPStatus::SUCCESS_EVENT;
                        this->message = "RungeKuttaFlow: Event triggered.";
                        this->success = true;
                        out.trim_chunk();
                        return out;
                    }
                    else
                    {
                        this->status = IVPStatus::SUCCESS_EVENT_BUT_TOL;
                        this->message = "RungeKuttaFlow: Event triggered, but integration tolerance was not met.";
                        this->success = true;
                        out.trim_chunk();
                        return out;
                    }
                }
            }

            const double next_error = solver.error_estimate;

            if (variable_time_step)
            {
                double dt_new_mult = 1.0;

                if (next_error < 1.0)
                {
                    // Step accepted

                    if (std::abs(next_error) < options.abstol)
                    {
                        // Prevent 0 to a power if error is 0.
                        dt_new_mult = options.large;
                    }
                    else
                    {
                        dt_new_mult = std::min(options.large, options.pessimist*std::pow(next_error, -1.0/(static_cast<double>(solver.order) + 1.0)));
                    }

                    dt_recommend = dt_new_mult*dt;
                    
                    // Obey min and max dt amounts
                    dt_recommend = std::min(dt_recommend, dt_max);
                    dt_recommend = std::max(dt_recommend, dt_min);
                    computestatus = 2;
                }
                else
                {
                    // Step rejected
                    if (dt - dt_min <= options.abstol)
                    {
                        // Accept step: dt is at lower bound. Continue solution as normal but warn.
                        this->tolerance_not_met = true;
                        computestatus = 2;
                    }
                    else
                    {
                        dt_new_mult = std::max(options.small, options.pessimist*std::pow(next_error, -1.0/(static_cast<double>(solver.order) + 1.0)));
                        dt = dt_new_mult*dt;
                        
                        // Obey min and max dt amounts
                        dt = std::min(dt, dt_max);
                        dt = std::max(dt, dt_min);
                        computestatus = 1;
                    }
                }
            }
            else
            {
                // Fixed time step
                computestatus = 2;
                if (!std::isnan(next_error) && next_error >= 1.0)
                {
                    // Accept step but signal warning (despite being fixed step).
                    this->tolerance_not_met = true;
                }
            }
        }
        else if (computestatus == 2)
        {
            // Step accepted. Advance solution.
            this->iterations += 1;
            t += dt;
            y.noalias() = ynext;
            out.add_data(t, y);
            
            if (!std::isnan(dt_save))
            {
                // If a saved dt is stored, reset it.
                dt = dt_save;
                dt_save = std::numeric_limits<double>::quiet_NaN();
            }
            else if (!std::isnan(dt_recommend))
            {
                // If a new dt was recommended, use it.
                dt = dt_recommend;
                dt_recommend = std::numeric_limits<double>::quiet_NaN();
            }

            // Advance tind if there are more in tspan
            if (t - tspan(tind) > -options.abstol)
            {
                tind += 1;
            }

            // Check all the termination conditions.
            if (t - tf > -options.abstol || tind >= tspan.size())
            {
                // Passed the last value in tspan
                if (!this->tolerance_not_met)
                {
                    this->status = IVPStatus::SUCCESS;
                    this->message = "RungeKuttaFlow: Final time reached.";
                    this->success = true;
                    out.trim_chunk();
                    return out;
                }
                else
                {
                    this->status = IVPStatus::SUCCESS_BUT_TOL;
                    this->message = "RungeKuttaFlow: Final time reached, but integration tolerance was not met.";
                    this->success = true;
                    out.trim_chunk();
                    return out;
                }
            }
            else if (options.max_iterations > 0 && this->iterations >= options.max_iterations)
            {
                // Max iterations exceeded.
                this->status = IVPStatus::ERROR_MAX_ITERATIONS;
                this->message = "RungeKuttaFlow: Max iterations reached.";
                this->success = false;
                out.trim_chunk();
                return out;
            }
            else
            {
                // Nothing was caught, keep running.
                if (t + dt - tspan(tind) >= -options.abstol)
                {
                    // Check if the next time step will cross the next time in tspan
                    // Save the current dt and set the next dt to end exactly on the next tspan value
                    dt_save = dt;
                    dt = tspan(tind) - t;
                }
                computestatus = 1;
            }
        }
    }

    // This should never be called
    this->status = IVPStatus::ERROR;
    this->message = "RungeKuttaFlow: Unknown error.";
    this->success = false;
    out.trim_chunk();
    return out;
}

}
