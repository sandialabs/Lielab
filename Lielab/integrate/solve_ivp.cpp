#include "solve_ivp.hpp"

#include "IVPMethods.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/testing.hpp"
#include "Lielab/utils.hpp"

#include <Eigen/Core>

#include <chrono>
#include <functional>
#include <limits>
#include <stdexcept>

namespace Lielab::integrate
{

IVPSolution solve_ivp(const EuclideanIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Eigen::VectorXd& y0, const IVPOptions options)
{
    /*!
    Solves a Euclidean IVP using the Runge-Kutta method.
    
    @param[in] dynamics A system defined on a Euclidean manifold.
    @param[in] tspan Independent variable.
    @param[in] y0 Initial condition.
    @param[in] options Additional options for the IVP solver.

    References
    ----------
    
    TODO
    ----
        - Preprocess values and error control
        - Flip tspan for descending values
    */

    using namespace Lielab::domain;
    IVPSolution out = IVPSolution();

    // Error checking
    lielab_assert(tspan.size() >= 2, "tspan must have size >= 2. Got " + std::to_string(tspan.size()) + ".");

    for (ptrdiff_t ii = 0; ii < tspan.size() - 1; ii++)
    {
        lielab_assert(tspan(ii+1) > tspan(ii), "tspan must be ascending.");
    }

    lielab_assert(options.dt_min < options.dt_max, "dt_min must be smaller than dt_max. Got " + std::to_string(options.dt_min) + " </= " + std::to_string(options.dt_max) + ".");

    const Eigen::VectorXd dy0 = dynamics.vectorfield(tspan(0), y0);

    lielab_assert(y0.size() == dy0.size(), "Dynamic system is not homogeneous. Expected: ... -> (R^" + std::to_string(y0.size()) + ") ->. Got: ... -> (R^" + std::to_string(dy0.size()) + ") ->.");

    // Preprocess options
    IVPOptions _options = options;
    if (_options.method == IVPMethod::Undefined)
    {
        _options.method = IVPMethod::RungeKutta;
    }

    // Do the computation
    if (_options.method == IVPMethod::RungeKutta || _options.method == IVPMethod::MuntheKaas)
    {
        // Main numerical method
        RungeKuttaFlow solver = RungeKuttaFlow();
        auto time_start = std::chrono::high_resolution_clock::now();
        out = solver(dynamics, tspan, y0, options);
        auto time_end = std::chrono::high_resolution_clock::now();

        // Postprocess
        const int n_eoms = static_cast<int>(out.ybar.cols());
        const int n_time = static_cast<int>(out.t.size());
        const CompositeManifold yzero = CompositeManifold({RN(n_eoms)});
        const CompositeAlgebra thetazero = CompositeAlgebra({rn(n_eoms)});

        if (!_options.rebase_every_step)
        {
            out.thetabar = out.ybar;
        }

        out.y = std::vector<CompositeManifold>(n_time, yzero);
        out.theta = std::vector<CompositeAlgebra>(n_time, thetazero);

        for (ptrdiff_t ii = 0; ii < n_time; ii++)
        {
            out.y[ii].unserialize(out.ybar.row(ii));
            out.theta[ii].set_vector(out.thetabar.row(ii));
        }

        out.time_to_solution = std::chrono::duration<double>(time_end - time_start).count();
        out.message = solver.message;
        out.success = solver.success;
        out.status = solver.status;
    }
    else if (_options.method == IVPMethod::CrouchGrossman)
    {
        throw std::runtime_error("CrouchGrossman method not implemented for Euclidean IVP Systems.");
    }
    else
    {
        throw std::runtime_error("Unknown method not implemented for Euclidean IVP Systems.");
    }

    return out;
}

IVPSolution solve_ivp(const HomogeneousIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Lielab::domain::CompositeManifold& y0, const IVPOptions options)
{
    /*!
    
    Solves a Homogeneous IVP using one of a few methods.

    CrouchGrossman
    --------------

    The Crouch-Grossman (CG) method [1].

    MuntheKaas (default)
    --------------------

    The Munthe-Kaas (MK) method [2-4].

    When `options.rebase_every_step = true`, this method will update the manifold position every timestep with
    theta = 0 as the original MK method intended. This yields a zero Lie algebra path for all time.

    When `options.rebase_every_step = false`, this method will start at theta = 0 for the initial condition, but
    it will let theta drift away from 0 resulting in a continuous, non-zero Lie algebra path.
    
    Parameters
    ----------

    @param[in] dynamics A system defined on a homogeneous manifold.
    @param[in] tspan Independent variable.
    @param[in] y0 Initial condition.
    @param[in] options Additional options for the IVP solver.

    References
    ----------
        [1] Crouch, Peter E., and R. Grossman. "Numerical integration of ordinary differential equations on manifolds."
                Journal of Nonlinear Science 3.1 (1993): 1-33.
        [2] Hans Munthe-Kaas. Lie-butcher theory for runge-kutta methods. BIT Numerical Mathematics, 35:572–587, 1995.
        [3] Hans Munthe-Kaas. Runge-kutta methods on lie groups. BIT Numerical Mathematics, 38:92–111, 1998.
        [4] Hans Munthe-Kaas. High order runge-kutta methods on manifolds. Applied Numerical Mathematics, 29(1):115–127, 1999
        [5] Hans Munthe-Kaas and Antonella Zanna. Numerical integration of differential equations on homogeneous manifolds.
                In Foundations of Computational Mathematics: Selected Papers of a Conference Held at
                Rio de Janeiro, January 1997, pages 305–315. Springer, 1997.

    */

    using namespace Lielab::domain;
    using Lielab::testing::check_topology;

    IVPSolution out = IVPSolution();

    // Error checking
    lielab_assert(tspan.size() >= 2, "tspan must have size >= 2. Got " + std::to_string(tspan.size()) + ".");

    for (ptrdiff_t ii = 0; ii < tspan.size() - 1; ii++)
    {
        lielab_assert(tspan(ii+1) > tspan(ii), "tspan must be ascending.");
    }

    lielab_assert(options.dt_min < options.dt_max, "dt_min must be smaller than dt_max. Got " + std::to_string(options.dt_min) + " </= " + std::to_string(options.dt_max) + ".");

    const CompositeAlgebra xi0 = dynamics.generator(tspan(0), y0);
    const CompositeAlgebra dy0 = dynamics.connection(0.0*xi0, xi0);
    const CompositeGroup Theta0 = dynamics.coordinates(dy0);
    const CompositeManifold ynext = dynamics.action(Theta0, y0);

    lielab_assert(check_topology(y0, ynext), "Dynamic system is not homogeneous. Expected: ... -> (" + y0.to_string() + "). Got: ... -> (" + ynext.to_string() + ").");

    // Preprocess options
    IVPOptions _options = options;
    if (_options.method == IVPMethod::Undefined)
    {
        _options.method = IVPMethod::MuntheKaas;
    }

    // Do the computation
    if (_options.method == IVPMethod::RungeKutta)
    {
        throw std::runtime_error("RungeKutta method not implemented for Homogeneous IVP Systems.");
    }
    else if (_options.method == IVPMethod::CrouchGrossman)
    {
        // Main numerical method
        CrouchGrossmanFlow solver;
        auto time_start = std::chrono::high_resolution_clock::now();
        out = solver(dynamics, tspan, y0, _options);
        auto time_end = std::chrono::high_resolution_clock::now();

        // Postprocessing
        const size_t n_t = out.t.size();
        out.y = std::vector<CompositeManifold>(n_t, y0);
        out.theta = std::vector<CompositeAlgebra>(n_t, 0.0*xi0);
        out.thetabar = Eigen::MatrixXd::Zero(n_t, xi0.get_vector().size());

        for (ptrdiff_t jj = 1; jj < out.t.size(); jj++)
        {
            out.y[jj].unserialize(out.ybar.row(jj));
        }

        out.time_to_solution = std::chrono::duration<double>(time_end - time_start).count();
        out.message = solver.message;
        out.success = solver.success;
        out.status = solver.status;
    }
    else if (_options.method == IVPMethod::MuntheKaas)
    {
        if (_options.rebase_every_step)
        {
            // Use the MK method as originally intended

            // Main numerical method
            MuntheKaasFlow solver;
            auto time_start = std::chrono::high_resolution_clock::now();
            out = solver(dynamics, tspan, y0, _options);
            auto time_end = std::chrono::high_resolution_clock::now();

            // Postprocessing
            const size_t n_t = out.t.size();
            out.y = std::vector<CompositeManifold>(n_t, y0);
            out.theta = std::vector<CompositeAlgebra>(n_t, 0.0*xi0);
            out.thetabar = Eigen::MatrixXd::Zero(n_t, xi0.get_vector().size());

            for (ptrdiff_t jj = 1; jj < out.t.size(); jj++)
            {
                out.y[jj].unserialize(out.ybar.row(jj));
            }

            out.time_to_solution = std::chrono::duration<double>(time_end - time_start).count();
            out.message = solver.message;
            out.success = solver.success;
            out.status = solver.status;
        }
        else
        {
            // Let theta drift away from 0

            // Preprocessing
            auto generator_wrapped = [&](const double t, const Eigen::VectorXd& thetabar)
            {
                CompositeAlgebra theta = 0.0*xi0; // Awkward statement forcing a copy.
                theta.set_vector(thetabar);
                const CompositeGroup Theta = dynamics.coordinates(theta);
                const CompositeManifold y = dynamics.action(Theta, y0);
                const CompositeAlgebra xi = dynamics.generator(t, y);
                const CompositeAlgebra dy = dynamics.connection(theta, xi);
                return dy.get_vector();
            };

            auto event_wrapped = [&](const double t, const Eigen::VectorXd& thetabar)
            {
                CompositeAlgebra theta = 0.0*xi0; // Awkward statement forcing a copy.
                theta.set_vector(thetabar);
                const CompositeGroup Theta = dynamics.coordinates(theta);
                const CompositeManifold y = dynamics.action(Theta, y0);
                return dynamics.event(t, y);
            };

            EuclideanIVPSystem MuntheKaasZannaDynamics(generator_wrapped);
            MuntheKaasZannaDynamics.event = event_wrapped;

            const Eigen::VectorXd thetabar0 = Eigen::VectorXd::Zero(xi0.get_dimension());

            IVPOptions rkoptions = _options;
            rkoptions.method = IVPMethod::RungeKutta;
            
            // Main numerical method
            out = solve_ivp(MuntheKaasZannaDynamics, tspan, thetabar0, rkoptions);

            // Postprocessing
            const size_t n_t = out.t.size();
            out.theta = std::vector<CompositeAlgebra>(n_t);
            out.y = std::vector<CompositeManifold>(n_t);
            out.ybar = Eigen::MatrixXd::Zero(n_t, y0.serialize().size());

            CompositeAlgebra thetaj = 0.0*xi0;
            thetaj.set_vector(out.thetabar.row(0));

            out.y[0] = y0;
            out.ybar.row(0) = y0.serialize();
            out.theta[0] = thetaj;

            for (ptrdiff_t jj = 1; jj < out.t.size(); jj++)
            {
                
                thetaj.set_vector(out.thetabar.row(jj));
                out.theta[jj] = 1.0*thetaj;
                const CompositeGroup Thetaj = dynamics.coordinates(thetaj);
                const CompositeManifold yj = dynamics.action(Thetaj, y0);
                out.y[jj] = yj;
                out.ybar.row(jj) = yj.serialize();
            }
        }
    }
    else
    {
        throw std::runtime_error("Unknown method for Homogeneous IVP Systems.");
    }

    return out;
}

}
