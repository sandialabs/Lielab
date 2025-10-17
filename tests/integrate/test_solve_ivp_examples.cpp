#include <functional>
#include <memory>
#include <numbers>

#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

/*
Tests solve_ivp can be used with "real" problems.
*/

TEST_CASE("solve_ivp_RNxRN_RN", "[integrate]")
{
    /*!
    Tests solve_ivp with the Lorenz equations.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;
    using Eigen::placeholders::last;

    const auto eoms = [](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const Eigen::VectorXd ybar = std::get<RN>(y.space[0]).serialize();
        
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

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0bar(0));
    CHECK(curve.ybar(0, 1) == y0bar(1));
    CHECK(curve.ybar(0, 2) == y0bar(2));

    CHECK(std::abs(curve.t(last) - 5.0) < TOL_FINE);
    CHECK(std::abs(curve.ybar(last, 0) - 15.230) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 1) + 0.797) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 2) + 1.473) < 1e-1);
}

TEST_CASE("solve_ivp_Composite_Coadjoint", "[integrate]")
{
    /*!
    Tests solve_ivp against a function with composite custom action (SE x GLR) x (SE x RN) -> (SE x RN).
    Also checks the integration of coadjoint actions.

    TODO
    ----
        - Rebasing every step will significantly improve accuracy
    */

    using namespace Lielab::domain;
    using namespace Lielab::functions;
    using namespace Lielab::integrate;
    using namespace Lielab::utils;
    using Eigen::placeholders::last;

    const auto eoms = [](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const double V = 1.0;
        const Eigen::VectorXd lambdabar = std::get<RN>(y.space[1]).serialize();
        const double u = -lambdabar(2);

        const se dx = se::from_vector({V, 0.0, u});
        const glr dlambdabar = -coad(dx);
        
        return CompositeAlgebra({dx, dlambdabar});
    };

    const auto action = [](const CompositeGroup& g, const CompositeManifold& y) -> CompositeManifold
    {
        const SE g0 = std::get<SE>(g.space[0]);
        const SE y0 = std::get<SE>(y.space[0]);
        const Eigen::MatrixXd coAdyhat = std::get<GLR>(g.space[1]).get_matrix();
        const Eigen::VectorXd lambdabar = std::get<RN>(y.space[1]).serialize();
        const RN lambdanext = RN::from_vector(coAdyhat*lambdabar);
        return CompositeManifold({y0*g0, lambdanext});
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

    Eigen::VectorXd tspan = linspace(0.0, 1.0, 100); // TODO:

    const CompositeManifold y0 = CompositeManifold({exp(x0), RN::from_vector(lambda0bar)});

    HomogeneousIVPSystem dynamics(eoms);
    dynamics.action = action;

    const ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 1) == -1.0);
    CHECK(curve.ybar(0, 3) == 1.0);
    CHECK(curve.ybar(0, 2) == x0bar(0));
    CHECK(curve.ybar(0, 5) == x0bar(1));
    CHECK(curve.ybar(0, 9) == lambda0bar(0));
    CHECK(curve.ybar(0, 10) == lambda0bar(1));
    CHECK(curve.ybar(0, 11) == lambda0bar(2));

    CHECK(std::abs(curve.t(last) - 1.0) < TOL_FINE);
    CHECK(std::abs(curve.ybar(last, 2) - 0.12732395447351627) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 5) + 0.0) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 9) - lambda0bar(0)) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 10) + lambda0bar(1)) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 11) - lambda0bar(2)) < 1e-1);
}
