#include "MuntheKaas.hpp"

#include "IVPCommon.hpp"

#include "Lielab/optimize.hpp"

#include <Eigen/Core>

#include <tuple>

namespace Lielab::integrate
{

MuntheKaas::MuntheKaas(const RungeKuttaCoefficients tableau)
{
    /*!
    * Instantiates a new MuntheKaas object with a specific RK method.
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

void MuntheKaas::estimate_error(const Eigen::VectorXd& yl, const Eigen::VectorXd& yh, const double dt)
{
    const Eigen::VectorXd scale = this->abstol + this->reltol*(yl.array().abs()).max(yh.array().abs());
    const int sz = static_cast<int>(yl.size());

    Eigen::VectorXd err0 = Eigen::VectorXd::Zero(sz);
    for (int ii = 0; ii < this->n; ii++)
    {
        err0 += this->K.row(ii)*this->e[ii]*dt;
    }

    const Eigen::VectorXd err = err0.array()/scale.array();
    this->error_estimate = err.norm()/std::sqrt(static_cast<double>(sz));
}

Lielab::domain::CompositeManifold MuntheKaas::operator()(const HomogeneousIVPSystem& dynamics, const Lielab::domain::CompositeManifold& y0, const double t0, const double dt)
{
    /*!
     * Primary usage of MuntheKaas.
     */

    const Lielab::domain::CompositeAlgebra dy0 = dynamics.generator(t0, y0);
    const Eigen::VectorXd dy0bar = dy0.get_vector();
    this->success = false;

    if (dy0bar.array().isNaN().any())
    {
        this->status = IVPStatus::ERROR_NANS_IN_VF;
        this->message = "MuntheKaas: NaNs in generator.";
        this->success = false;
        return y0;
    }

    if (dy0bar.array().isInf().any())
    {
        this->status = IVPStatus::ERROR_INFS_IN_VF;
        this->message = "MuntheKaas: Infs in generator.";
        this->success = false;
        return y0;
    }

    const int n_eoms = static_cast<int>(dy0bar.size());

    if (n_eoms == 0) return Lielab::domain::CompositeManifold();

    this->K = Eigen::MatrixXd::Zero(this->n, n_eoms);

    double next_t = t0;
    Lielab::domain::CompositeAlgebra next_theta = 0.0*dy0;

    int max_iter = 50;

    if (!this->implicit)
    {
        for (int stage = 0; stage < this->n; stage++)
        {
            if (stage == 0)
            {
                next_t = t0;
                next_theta = 0.0*dy0;
            }
            else
            {
                next_t = t0 + dt*this->C(stage);
                Eigen::VectorXd _next_theta = Eigen::VectorXd::Zero(n_eoms);
                for (int jj = 0; jj < stage; jj++)
                {
                    _next_theta += dt*this->A(stage, jj)*this->K.row(jj);
                }
                next_theta.set_vector(_next_theta);
            }

            this->K.row(stage) = dynamics.connection(next_theta, dynamics.generator(next_t, dynamics.action(dynamics.coordinates(next_theta), y0))).get_vector();
        }

        if (this->K.array().isNaN().any())
        {
            this->status = IVPStatus::ERROR_NANS_IN_VF;
            this->message = "MuntheKaas: NaNs in generator.";
            this->success = false;
            return y0;
        }

        if (this->K.array().isInf().any())
        {
            this->status = IVPStatus::ERROR_INFS_IN_VF;
            this->message = "MuntheKaas: Infs in generator.";
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
                Eigen::VectorXd _next_theta = Eigen::VectorXd::Zero(n_eoms);
                for (int jj = 0; jj < this->n; jj++)
                {
                    _next_theta += dt*this->A(stage, jj)*this->K.row(jj);
                }
                next_theta.set_vector(_next_theta);

                Knew.row(stage) = dynamics.connection(next_theta, dynamics.generator(next_t, dynamics.action(dynamics.coordinates(next_theta), y0))).get_vector();
            }

            err = std::sqrt((Knew - this->K).array().square().sum());
            this->K = Knew;

            if (this->K.array().isNaN().any())
            {
                this->status = IVPStatus::ERROR_NANS_IN_VF;
                this->message = "MuntheKaas: NaNs in generator.";
                this->success = false;
                return y0;
            }

            if (this->K.array().isInf().any())
            {
                this->status = IVPStatus::ERROR_INFS_IN_VF;
                this->message = "MuntheKaas: Infs in generator.";
                this->success = false;
                return y0;
            }

            iteration++;
        }
        
        // TODO: Error if max iter here
    }


    next_theta *= 0.0;
    Lielab::domain::CompositeAlgebra next_theta2 = 0.0*next_theta;
    Eigen::VectorXd _next_theta = Eigen::VectorXd::Zero(n_eoms);
    for (int ii = 0; ii < this->n; ii++)
    {
        _next_theta += dt*this->B(ii)*this->K.row(ii);
    }
    next_theta.set_vector(_next_theta);

    // Postprocess
    this->error_estimate = 0.0;
    if (this->can_variable_step == true)
    {
        Eigen::VectorXd _next_theta2 = Eigen::VectorXd::Zero(n_eoms);
        for (int ii = 0; ii < this->n; ii++)
        {
            _next_theta2 += dt*this->Bhat(ii)*this->K.row(ii);
        }
        this->estimate_error(_next_theta, _next_theta2, dt);
    }

    this->status = IVPStatus::SUCCESS;
    this->success = true;
    this->message = "MuntheKaas: Terminated successfully.";

    return dynamics.action(dynamics.coordinates(next_theta), y0);
}

/*
* The Flow class. Numerically calculates the integral curve of a vector field through y0 from t0 to tf.
*/

IVPSolution MuntheKaasFlow::operator()(const HomogeneousIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Lielab::domain::CompositeManifold& y0, const IVPOptions& options)
{
    /*!
    * The main evaluation method for the Flow class.
    */

    using namespace Lielab::optimize;

    // Preprocessing
    const ptrdiff_t n_states = y0.serialize().size();

    IVPSolution out(n_states);
    this->iterations = 0;
    this->message = "";
    this->success = false;
    this->status = IVPStatus::RUNNING;

    MuntheKaas solver = MuntheKaas(options.coefficients);

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
    Lielab::domain::CompositeManifold y = y0;
    Lielab::domain::CompositeManifold ynext = y0;

    out.add_data(tspan(0), y.serialize());

    // Main loop. Run until algorithm says it's done
    while (computestatus == 1 || computestatus == 2)
    {
        if (computestatus == 1)
        {
            // Construct next possible step
            ynext = solver.operator()(dynamics, y, t, dt);
            
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
                        const Lielab::domain::CompositeManifold _ynext = solver.operator()(dynamics, y, t, x(0));
                        return std::abs(dynamics.event(t + x(0), _ynext));
                    };

                    EuclideanExtremizationSystem system(fun);
                    system.lower_bound = Lielab::utils::to_VectorXd({std::numeric_limits<double>::epsilon()});
                    system.upper_bound = Lielab::utils::to_VectorXd({dt});
                    const Eigen::VectorXd xopt = search(system, optimoptions);
                    dt = xopt(0);

                    ynext = solver.operator()(dynamics, y, t, dt);
                    
                    this->iterations += 1;
                    t += dt;
                    y = ynext;
                    out.add_data(t, y.serialize());

                    // Exit
                    if (!this->tolerance_not_met)
                    {
                        this->status = IVPStatus::SUCCESS_EVENT;
                        this->message = "MuntheKaasFlow: Event triggered.";
                        this->success = true;
                        out.trim_chunk();
                        return out;
                    }
                    else
                    {
                        this->status = IVPStatus::SUCCESS_EVENT_BUT_TOL;
                        this->message = "MuntheKaasFlow: Event triggered, but integration tolerance was not met.";
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
            y = ynext;
            out.add_data(t, y.serialize());
            
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
                    this->message = "MuntheKaasFlow: Final time reached.";
                    this->success = true;
                    out.trim_chunk();
                    return out;
                }
                else
                {
                    this->status = IVPStatus::SUCCESS_BUT_TOL;
                    this->message = "MuntheKaasFlow: Final time reached, but integration tolerance was not met.";
                    this->success = true;
                    out.trim_chunk();
                    return out;
                }
            }
            else if (options.max_iterations > 0 && this->iterations >= options.max_iterations)
            {
                // Max iterations exceeded.
                this->status = IVPStatus::ERROR_MAX_ITERATIONS;
                this->message = "MuntheKaasFlow: Max iterations reached.";
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

    // This should never get executed
    this->status = IVPStatus::ERROR;
    this->message = "MuntheKaasFlow: Unknown error.";
    this->success = false;
    out.trim_chunk();
    return out;
}

}
