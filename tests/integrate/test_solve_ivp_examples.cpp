#include <functional>
#include <memory>
#include <numbers>

#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

/*
Tests solve_ivp can be used with "real" problems.
*/

TEST_CASE("solve_ivp_Lorenz_RNxRN_RN", "[integrate]")
{
    /*!
    Tests solve_ivp with the Lorenz equations.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;
    using Eigen::placeholders::last;

    const auto eoms = [](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        
        Eigen::VectorXd dybar(3);
        const double sigma = 10.0;
        const double rho = 28.0;
        const double b1 = 8.0;
        const double b2 = 3.0;
        const double beta = b1/b2;

        dybar(0) = -beta*ybar(0) + ybar(1)*ybar(2);
        dybar(1) = -sigma*ybar(1) + sigma*ybar(2);
        dybar(2) = -ybar(0)*ybar(1) + rho*ybar(1) - ybar(2);
        return CompositeAlgebra({rn::from_vector(dybar)});
    };

    Eigen::VectorXd y0bar(3);
    y0bar(0) = 25.0;
    y0bar(1) = 0.0;
    y0bar(2) = -20.0;
    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = 5.0;

    HomogeneousIVPSystem dynamics(eoms);

    IVPSolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0bar(0));
    CHECK(curve.ybar(0, 1) == y0bar(1));
    CHECK(curve.ybar(0, 2) == y0bar(2));

    CHECK(std::abs(curve.t(last) - 5.0) < TOL_FINE);
    CHECK(std::abs(curve.ybar(last, 0) - 15.230) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 1) + 0.797) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 2) + 1.473) < 1e-1);
}

TEST_CASE("solve_ivp_Pleiades", "[integrate]")
{
    /*
    N-body problem with n=7 "Pleiades Problem".

    This could probably get higher accuracy with GLR actions but needs rebasing every step in MK.

    References
    ----------
        [1] Hairer, Ernst, Gerhard Wanner, and Syvert P. Nørsett. Solving ordinary differential equations I:
            Nonstiff problems. Berlin, Heidelberg: Springer Berlin Heidelberg, 1993.
        [2] http://archimede.dm.uniba.it/~testset/report/plei.pdf
    */
    
    using namespace Lielab::integrate;
    using namespace Lielab::utils;
    using Eigen::placeholders::last;

    Eigen::VectorXd m = to_VectorXd({1, 2, 3, 4, 5, 6, 7});

    auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
    {
        const int n = static_cast<int>(y.size()/4);

        const Eigen::VectorXd px = y(Eigen::seqN(0, n));
        const Eigen::VectorXd py = y(Eigen::seqN(n, n));
        const Eigen::VectorXd vx = y(Eigen::seqN(2*n, n));
        const Eigen::VectorXd vy = y(Eigen::seqN(3*n, n));

        Eigen::VectorXd pxd = vx;
        Eigen::VectorXd pyd = vy;
        Eigen::VectorXd vxd = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd vyd = Eigen::VectorXd::Zero(n);

        for (ptrdiff_t ii = 0; ii < n; ii++)
        {
            for (ptrdiff_t jj = 0; jj < n; jj++)
            {
                if (ii != jj)
                {
                    const double rad = std::pow(std::pow(px[ii] - px[jj], 2.0) + std::pow(py[ii] - py[jj], 2.0), 1.5);
                    vxd(ii) += m[jj]*(px[jj] - px[ii])/rad;
                    vyd(ii) += m[jj]*(py[jj] - py[ii])/rad;
                }
            }
        }
        
        return concatenate({pxd, pyd, vxd, vyd});
    };
    
    EuclideanIVPSystem dynamics(eoms);
    IVPOptions options = IVPOptions();
    options.abstol = 1e-10;
    options.reltol = 1e-10;

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = 3.0;

    const Eigen::VectorXd y0 = to_VectorXd({3.0, 3.0, -1.0, -3.0, 2.0, -2.0, 2.0, 3.0, -3.0, 2.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.75, -1.5, 0.0, 0.0, 0.0, -1.25, 1.0, 0.0, 0.0});
    IVPSolution path = solve_ivp(dynamics, tspan, y0, options);

    const double tol = 1e-7;

    CHECK(std::abs(path.ybar(last, 0) - 0.3706139143970502) < tol);
    CHECK(std::abs(path.ybar(last, 1) - 3.237284092057233) < tol);
    CHECK(std::abs(path.ybar(last, 2) + 3.222559032418324) < tol);
    CHECK(std::abs(path.ybar(last, 3) - 0.6597091455775310) < tol);
    CHECK(std::abs(path.ybar(last, 4) - 0.3425581707156584) < tol);
    CHECK(std::abs(path.ybar(last, 5) - 1.562172101400631) < tol);
    CHECK(std::abs(path.ybar(last, 6) + 0.7003092922212495) < tol);

    CHECK(std::abs(path.ybar(last, 7) + 3.943437585517392) < tol);
    CHECK(std::abs(path.ybar(last, 8) + 3.271380973972550) < tol);
    CHECK(std::abs(path.ybar(last, 9) - 5.225081843456543) < tol);
    CHECK(std::abs(path.ybar(last, 10) + 2.590612434977470) < tol);
    CHECK(std::abs(path.ybar(last, 11) - 1.198213693392275) < tol);
    CHECK(std::abs(path.ybar(last, 12) + 0.2429682344935824) < tol);
    CHECK(std::abs(path.ybar(last, 13) - 1.091449240428980) < tol);

    CHECK(std::abs(path.ybar(last, 14) - 3.417003806314313) < tol);
    CHECK(std::abs(path.ybar(last, 15) - 1.354584501625501) < tol);
    CHECK(std::abs(path.ybar(last, 16) + 2.590065597810775) < tol);
    CHECK(std::abs(path.ybar(last, 17) - 2.025053734714242) < tol);
    CHECK(std::abs(path.ybar(last, 18) + 1.155815100160448) < tol);
    CHECK(std::abs(path.ybar(last, 19) + 0.8072988170223021) < tol);
    CHECK(std::abs(path.ybar(last, 20) - 0.5952396354208710) < tol);

    CHECK(std::abs(path.ybar(last, 21) + 3.741244961234010) < tol);
    CHECK(std::abs(path.ybar(last, 22) - 0.3773459685750630) < tol);
    CHECK(std::abs(path.ybar(last, 23) - 0.9386858869551073) < tol);
    CHECK(std::abs(path.ybar(last, 24) - 0.3667922227200571) < tol);
    CHECK(std::abs(path.ybar(last, 25) + 0.3474046353808490) < tol);
    CHECK(std::abs(path.ybar(last, 26) - 2.344915448180937) < tol);
    CHECK(std::abs(path.ybar(last, 27) + 1.947020434263292) < tol);
}

TEST_CASE("solve_ivp_Particle_Magnetic", "[integrate]")
{
    /*
    Solves the motion of a charged particle in a magnetic dipole field. Ex 11.3 from Ref 1. Note that
    the authors construct the action as (GL(6), R^6) -> R^6 and do exponential coordinates on
    the entire gl(6)/GL(6). Numerically this isn't very accurate, especially as theta drifts away from the
    0 of the algebra. So I rewrote the action as (R^3 x SO(3), R^3 x R^3) -> (R^3 x R^3) since SO(3)
    is much more well behaved. Resulting path is a particle trapped in the van Allen belt.

    References
    ----------
        [1] Iserles, Arieh, Hans Z. Munthe-Kaas, Syvert P. Nørsett, and Antonella Zanna. "Lie-group methods." Acta numerica 9 (2000): 215-365.
    */
    
    using namespace Lielab::domain;
    using namespace Lielab::integrate;
    using namespace Lielab::utils;
    using Eigen::placeholders::last;

    const double pi = std::numbers::pi;

    const Eigen::VectorXd tspan = linspace<double>(0.0, 500.0, 2);
    const CompositeManifold y0 = CompositeManifold({RN::from_vector({0.0, -2.5, 0.0}), RN::from_vector({0.0, 0.0, 0.012})});

    auto vectorfield = [pi](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const RN& y1 = y[1];
        const Eigen::VectorXd p = y0.serialize();
        const Eigen::VectorXd v = y1.serialize();
        
        const double p02 = std::pow(p(0), 2.0);
        const double p12 = std::pow(p(1), 2.0);
        const double p22 = std::pow(p(2), 2.0);

        const double azimuth = std::atan2(p(1), p(0));
        const double elevation = std::atan2(p(2), std::sqrt(p02 + p12));
        const double rad = std::sqrt(p02 + p12 + p22);

        // m = [0, 0, 1] // Magnetic dipole vector
        const double rad3 = std::pow(rad, 3.0);
        const double br = 2.0*std::sin(elevation)/rad3;
        const double bx = br*std::cos(elevation)*std::cos(azimuth);
        const double by = br*std::cos(elevation)*std::sin(azimuth);
        const double bz = br*std::sin(elevation);
        const double btheta = std::cos(elevation)/rad3;
        const double bthx = btheta*std::cos(elevation - pi/2.0)*std::cos(azimuth);
        const double bthy = btheta*std::cos(elevation - pi/2.0)*std::sin(azimuth);
        const double bthz = btheta*std::sin(elevation - pi/2.0);
        const double b1 = bx + bthx;
        const double b2 = by + bthy;
        const double b3 = bz + bthz;

        return CompositeAlgebra({rn::from_vector(v), so::from_vector({b1,b2,b3})});
    };
    
    auto action = [](const CompositeGroup& g, const CompositeManifold& y) -> CompositeManifold
    {
        const RN& g0 = g[0];
        const SO& g1 = g[1];
        const RN& y0 = y[0];
        const RN& y1 = y[1];

        const Eigen::VectorXd pnext = g0.serialize() + y0.serialize();
        const Eigen::VectorXd vnext = g1.get_matrix()*y1.serialize();
        return CompositeManifold({RN::from_vector(pnext), RN::from_vector(vnext)});
    };
    
    HomogeneousIVPSystem dynamics(vectorfield);
    dynamics.action = action;

    IVPOptions options;
    options.dt_max = 1.0;

    IVPSolution path = solve_ivp(dynamics, tspan, y0, options);

    REQUIRE(path.success == true);

    CHECK(std::abs(path.ybar(last, 0) - -0.527505080797981) < 1e-3);
    CHECK(std::abs(path.ybar(last, 1) - -2.434436925360713) < 1e-3);
    CHECK(std::abs(path.ybar(last, 2) - -0.145842164780228) < 1e-3);
    CHECK(std::abs(path.ybar(last, 3) - 0.000006752678154) < 1e-3);
    CHECK(std::abs(path.ybar(last, 4) - 0.001167463126786) < 1e-3);
    CHECK(std::abs(path.ybar(last, 5) - -0.011943072646892) < 1e-3);

    // MK integration perfectly preserves Hamiltonian = |v|^2 since acceleration is orthogonal to velocity.
    const RN& v0 = path.y[0][1];
    const RN& vf = path.y.back()[1];
    const double h0 = v0.serialize().norm();
    const double hf = vf.serialize().norm();
    CHECK(std::abs(h0 - hf) < 1e-16);
}

TEST_CASE("solve_ivp_Composite_Coadjoint", "[integrate]")
{
    /*!
    Tests solve_ivp against a function with composite custom action (SE x GLR) x (SE x RN) -> (SE x RN).
    Also checks the integration of coadjoint actions and custom actions.
    */

    using namespace Lielab::domain;
    using namespace Lielab::functions;
    using namespace Lielab::integrate;
    using namespace Lielab::utils;
    using Eigen::placeholders::last;

    const auto eoms = [](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& lambda = y[1];
        const double V = 1.0;
        const Eigen::VectorXd lambdabar = lambda.serialize();
        const double u = -lambdabar(2);

        const se dx = se::from_vector({V, 0.0, u});
        const glr dlambdabar = -coad(dx);
        
        return CompositeAlgebra({dx, dlambdabar});
    };

    const auto action = [](const CompositeGroup& g, const CompositeManifold& y) -> CompositeManifold
    {
        const SE& g0 = g[0];
        const GLR& g1 = g[1];
        const SE& y0 = y[0];
        const RN& y1 = y[1];
        const Eigen::MatrixXd coAdyhat = g1.get_matrix();
        const Eigen::VectorXd lambdabar = y1.serialize();
        const RN lambdanext = RN::from_vector(coAdyhat*lambdabar);
        return CompositeManifold({y0*g0, lambdanext});
    };

    const auto connection = [](const CompositeAlgebra& theta, const CompositeAlgebra& xi) -> CompositeAlgebra
    {
        const se& theta0 = theta[0];
        const se& xi0 = xi[0];
        const glr& theta1 = theta[1];
        const glr& xi1 = xi[1];
        return CompositeAlgebra({dexpinv(-theta0, xi0), dexpinv(theta1, xi1)});
    };

    Eigen::VectorXd x0bar(3);
    x0bar(0) = 0.0;
    x0bar(1) = 0.0;
    x0bar(2) = std::numbers::pi_v<double>/2.0;
    const se x0 = se::from_vector(x0bar);

    Eigen::VectorXd lambda0bar(3);
    lambda0bar(0) = 1.15407533e-03;
    lambda0bar(1) = -3.17495766e+01;
    lambda0bar(2) = -4.41935411e+00;

    Eigen::VectorXd tspan = linspace(0.0, 1.0, 2);

    const CompositeManifold y0 = CompositeManifold({exp(x0), RN::from_vector(lambda0bar)});

    HomogeneousIVPSystem dynamics(eoms);
    dynamics.action = action;
    dynamics.connection = connection;

    IVPOptions options;
    options.dt_max = 0.1;

    const IVPSolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == x0bar(0));
    CHECK(curve.ybar(0, 1) == x0bar(1));
    CHECK(curve.ybar(0, 3) == -1.0);
    CHECK(curve.ybar(0, 4) == 1.0);
    CHECK(curve.ybar(0, 6) == lambda0bar(0));
    CHECK(curve.ybar(0, 7) == lambda0bar(1));
    CHECK(curve.ybar(0, 8) == lambda0bar(2));

    CHECK(std::abs(curve.t(last) - 1.0) < TOL_FINE);
    CHECK(std::abs(curve.ybar(last, 0) - 0.12732395447351627) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 1) + 0.0) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 6) - lambda0bar(0)) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 7) + lambda0bar(1)) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 8) - lambda0bar(2)) < 1e-1);
}

TEST_CASE("solve_ivp_AEIF", "[integrate]")
{
    /*!
    Tests solve_ivp with the Adaptive Exponential Integrate and Fire model. Primarily checks
    robustness of event handling.

    References
    ----------
        [1] https://en.wikipedia.org/wiki/Exponential_integrate-and-fire
        [2] https://nest-simulator.readthedocs.io/en/stable/model_details/aeif_models_implementation.html
    */

    using namespace Lielab::integrate;
    using namespace Lielab::utils;

    const double V_reset = -58.0;
    const double V_peak = 0.0;
    const double V_th = -50.0;
    const double I_e = 420.0;
    const double g_L = 11.0;
    const double tau_w = 300.0;
    const double E_L = -70.0;
    const double Delta_T = 2.0;
    const double a = 3.0;
    const double C_m = 200.0;

    auto eoms = [&](const double t, const Eigen::VectorXd& y)
    {
        const double v = std::min(y(0), V_peak);
        const double w = y(1);

        double Ispike = 0.0;

        if (Delta_T != 0.0)
        {
            Ispike = g_L * Delta_T * std::exp((v - V_th) / Delta_T);
        }

        const double dv = (-g_L * (v - E_L) + Ispike - w + I_e) / C_m;
        const double dw = (a * (v - E_L) - w) / tau_w;

        return to_VectorXd({dv, dw});
    };

    auto event = [&](const double t, const Eigen::VectorXd& y)
    {
        return V_peak - y(0);
    };

    EuclideanIVPSystem dynamics(eoms);
    dynamics.event = event;

    IVPOptions options;
    options.dt_max = 0.01;

    const double tf = 100.0;
    bool first_run = true;
    bool tf_reached = false;
    int n_events = 0;

    IVPSolution sol;
    Eigen::VectorXd tspan;
    Eigen::VectorXd y0;

    while (!tf_reached)
    {
        if (first_run)
        {
            tspan = to_VectorXd({0.0, tf});
            y0 = to_VectorXd({E_L, 5.0});
        }
        else
        {
            const ptrdiff_t n_t = sol.t.size();
            tspan = to_VectorXd({sol.t(n_t - 1), tf});
            y0 = to_VectorXd({V_reset, sol.ybar(n_t - 1, 1)});
        }
        
        const IVPSolution seg = solve_ivp(dynamics, tspan, y0, options);

        if (seg.status == IVPStatus::SUCCESS_EVENT || seg.status == IVPStatus::SUCCESS_EVENT_BUT_TOL)
        {
            const ptrdiff_t n_t = seg.t.size();
            CHECK(std::abs(seg.ybar(n_t - 1, 0) - V_peak) < 0.1);
            n_events++;
        }
        
        if (first_run)
        {
            sol = seg;
        }
        else
        {
            sol.t = concatenate({sol.t, seg.t});
            sol.ybar = vertical_stack({sol.ybar, seg.ybar});
            sol.thetabar = vertical_stack({sol.thetabar, seg.thetabar});
            sol.y.insert(sol.y.end(), seg.y.begin(), seg.y.end());
            sol.theta.insert(sol.theta.end(), seg.theta.begin(), seg.theta.end());
        }

        const ptrdiff_t n_total = sol.t.size();
        if (sol.t(n_total - 1) >= tf)
        {
            tf_reached = true;
        }

        first_run = false;
    }

    CHECK(n_events == 7);
    const ptrdiff_t n_total = sol.t.size();
    CHECK(std::abs(sol.t(n_total - 1) - tf) < 1e-15);
}
