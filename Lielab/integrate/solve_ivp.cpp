#include "solve_ivp.hpp"

#include "Coefficients.hpp"
#include "ODESolution.hpp"

#include "IVPMethods/IVPSettings.hpp"
#include "IVPMethods/RungeKutta.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/utils.hpp"

#include <functional>
#include <limits>

namespace Lielab::integrate
{

ODESolution solve_ivp(const EuclideanIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Eigen::VectorXd& y0, const IVPOptions options)
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

    */

    RungeKuttaFlow F = RungeKuttaFlow();

    F.dt = options.dt_initial;
    F.dt_min = options.dt_min;
    F.dt_max = options.dt_max;

    F.reltol = options.reltol;
    F.abstol = options.abstol;

    F.method->reltol = options.reltol;
    F.method->abstol = options.abstol;

    F.variable_time_step = options.variable_time_step;

    ODESolution out = F(dynamics, tspan, y0, options);

    if (out.status == 0)
    {
        out.message = "Method converged.";
        out.success = true;
    }
    else if (out.status == -1)
    {
        out.message = "NaNs in vectorfield.";
        out.success = false;
    }
    else if (out.status == -2)
    {
        out.message = "Infs in vectorfield.";
        out.success = false;
    }

    return out;
}

ODESolution solve_ivp(const HomogeneousIVPSystem& dynamics, const Eigen::VectorXd& tspan, const Lielab::domain::CompositeManifold& y0, const IVPOptions options)
{
    /*!
    Solves a Homogeneous IVP using the Munthe-Kaas method [1-3]. Re-uses the Euclidean Runge-Kutta solver [4].
    
    @param[in] dynamics A system defined on a homogeneous manifold.
    @param[in] tspan Independent variable.
    @param[in] y0 Initial condition.
    @param[in] options Additional options for the IVP solver.

    References
    ----------
        [1] Hans Munthe-Kaas. Lie-butcher theory for runge-kutta methods. BIT Numerical Mathematics, 35:572–587, 1995.
        [2] Hans Munthe-Kaas. Runge-kutta methods on lie groups. BIT Numerical Mathematics, 38:92–111, 1998.
        [3] Hans Munthe-Kaas. High order runge-kutta methods on manifolds. Applied Numerical Mathematics, 29(1):115–127, 1999
        [4] Hans Munthe-Kaas and Antonella Zanna. Numerical integration of differential equations on homogeneous manifolds.
                In Foundations of Computational Mathematics: Selected Papers of a Conference Held at
                Rio de Janeiro, January 1997, pages 305–315. Springer, 1997.
    
    TODO
    ----
        - Reusing the RK solver tends to cause theta to drift away from the 0 of the algebra giving certain coordinate
          mappings numerical difficulty. Rebasing theta about 0 every time step will mitigate this.

    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;

    // TODO: Error check here

    const CompositeAlgebra xi0 = dynamics.vectorfield(tspan(0), y0);

    std::vector<ODESolution> rksols;
    std::vector<CompositeManifold> y0seg;
    y0seg.push_back(y0);

    for (size_t ii = 0; ii < tspan.size() - 1; ii++)
    {
        auto vf_wrapped = [&](const double t, const Eigen::VectorXd& thetabar)
        {
            CompositeAlgebra theta = 0.0*xi0; // Awkward statement forcing a copy.
            theta.set_vector(thetabar);
            const CompositeGroup Theta = dynamics.coordinates(theta);
            const CompositeManifold y = dynamics.action(Theta, y0seg[ii]);
            const CompositeAlgebra xi = dynamics.vectorfield(t, y);
            const CompositeAlgebra dy = dynamics.connection(theta, xi);
            return dy.get_vector();
        };

        auto event_wrapped = [&](const double t, const Eigen::VectorXd& thetabar)
        {
            CompositeAlgebra theta = 0.0*xi0; // Awkward statement forcing a copy.
            theta.set_vector(thetabar);
            const CompositeGroup Theta = dynamics.coordinates(theta);
            const CompositeManifold y = dynamics.action(Theta, y0seg[ii]);
            return dynamics.event(t, y);
        };

        Eigen::VectorXd tseg(2);
        tseg(0) = tspan(ii);
        tseg(1) = tspan(ii+1);

        EuclideanIVPSystem MuntheKaasZannaDynamics;
        MuntheKaasZannaDynamics.vectorfield = vf_wrapped;
        MuntheKaasZannaDynamics.event = event_wrapped;

        const Eigen::VectorXd thetabar0 = Eigen::VectorXd::Zero(xi0.get_dimension());

        ODESolution segment = solve_ivp(MuntheKaasZannaDynamics, tseg, thetabar0, options);
        const size_t n_t = segment.t.size();
        segment.thetabar = 1.0*segment.ybar;
        segment.theta = std::vector<CompositeAlgebra>(n_t);
        segment.y = std::vector<CompositeManifold>(n_t);
        segment.ybar = Eigen::MatrixXd::Zero(n_t, y0.serialize().size());

        for (size_t jj = 0; jj < segment.t.size(); jj++)
        {
            CompositeAlgebra thetaj = 0.0*xi0;
            thetaj.set_vector(segment.thetabar.row(jj));
            segment.theta[jj] = 1.0*thetaj;
            const CompositeGroup Thetaj = dynamics.coordinates(thetaj);
            const CompositeManifold yj = dynamics.action(Thetaj, y0seg[ii]);
            segment.y[jj] = yj;
            segment.ybar.row(jj) = yj.serialize();
        }

        rksols.push_back(segment);
        y0seg.push_back(segment.y.back());
    }

    ODESolution out = rksols[0];

    for (size_t ii = 0; ii < tspan.size() - 2; ii++)
    {
        const size_t n_t = rksols[ii+1].t.size();
        const Eigen::VectorXd next_t = rksols[ii+1].t(Eigen::seqN(1, n_t - 1));
        const Eigen::MatrixXd next_thetabar = rksols[ii+1].thetabar(Eigen::seqN(1, n_t - 1), Eigen::all);
        const Eigen::MatrixXd next_ybar = rksols[ii+1].ybar(Eigen::seqN(1, n_t - 1), Eigen::all);
        const std::vector<CompositeAlgebra> next_theta(rksols[ii+1].theta.begin() + 1, rksols[ii+1].theta.end());
        const std::vector<CompositeManifold> next_y(rksols[ii+1].y.begin() + 1, rksols[ii+1].y.end());

        out.t = Lielab::utils::concatenate({out.t, next_t});
        out.thetabar = Lielab::utils::vertical_stack({out.thetabar, next_thetabar});
        out.ybar = Lielab::utils::vertical_stack({out.ybar, next_ybar});
        out.theta.insert(out.theta.end(), next_theta.begin(), next_theta.end());
        out.y.insert(out.y.end(), next_y.begin(), next_y.end());
    }

    return out;
}

}
