#include <functional>
#include <memory>
#include <numbers>

#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

// Exponential decay tests

TEST_CASE("solve_ivp_1_euclidean", "[integrate]")
{
    /*!
    * Tests solve_ivp against a classical problem with known solution.
    */

    using namespace Lielab::integrate;

    const double m = -0.5;
    const double tf = 10.0;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    // Do not modify anything below this line.

    const auto eoms = [m](const double t, const Eigen::VectorXd& y)
    {
        Eigen::VectorXd dy(4);

        dy(0) = m*y(0);
        dy(1) = m*y(1);
        dy(2) = m*y(2);
        dy(3) = m*y(3);
        return dy;
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    EuclideanIVPSystem dynamics;
    dynamics.vectorfield = eoms;

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0(0));
    CHECK(curve.ybar(0, 1) == y0(1));
    CHECK(curve.ybar(0, 2) == y0(2));
    CHECK(curve.ybar(0, 3) == y0(3));

    const size_t nrows = curve.ybar.rows();

    CHECK(std::abs(curve.t(nrows-1) - tf) < TOL_FINE);
    CHECK(std::abs(curve.ybar(nrows-1, 0) - std::exp(m*tf)*y0(0)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 1) - std::exp(m*tf)*y0(1)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 2) - std::exp(m*tf)*y0(2)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 3) - std::exp(m*tf)*y0(3)) < 1e-5);
}

TEST_CASE("solve_ivp_1_euclidean_event", "[integrate]")
{
    /*!
    * Tests solve_ivp against a classical problem with known solution.
    */

    using namespace Lielab::integrate;

    const double m = -0.5;
    const double tf = 10.0;
    const double y2cross = 3.5;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    // Do not modify anything below this line.

    const auto eoms = [m](const double t, const Eigen::VectorXd& y)
    {
        Eigen::VectorXd dy(4);

        dy(0) = m*y(0);
        dy(1) = m*y(1);
        dy(2) = m*y(2);
        dy(3) = m*y(3);
        return dy;
    };

    const auto event = [y2cross](const double t, const Eigen::VectorXd& y)
    {
        return y(2) - y2cross;
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    EuclideanIVPSystem dynamics;
    dynamics.vectorfield = eoms;
    dynamics.event = event;

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0(0));
    CHECK(curve.ybar(0, 1) == y0(1));
    CHECK(curve.ybar(0, 2) == y0(2));
    CHECK(curve.ybar(0, 3) == y0(3));

    const size_t nrows = curve.ybar.rows();
    const double tcross = std::log(y2cross/y0(2))/m;

    CHECK(std::abs(curve.t(nrows-1) - tcross) < TOL_FINE);
    CHECK(std::abs(curve.ybar(nrows-1, 0) - std::exp(m*tcross)*y0(0)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 1) - std::exp(m*tcross)*y0(1)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 2) - std::exp(m*tcross)*y0(2)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 3) - std::exp(m*tcross)*y0(3)) < 1e-5);
}

TEST_CASE("solve_ivp_1", "[integrate]")
{
    /*!
    * Tests solve_ivp against a classical problem with known solution.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;

    const double m = -0.5;
    const double tf = 10.0;

    Eigen::VectorXd y0bar(4);
    y0bar(0) = 2.0;
    y0bar(1) = 4.0;
    y0bar(2) = 6.0;
    y0bar(3) = -8.0;

    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    // Do not modify anything below this line.

    const auto eoms = [m](const double t, const CompositeManifold& y)
    {
        const Eigen::VectorXd ybar = std::get<RN>(y.space[0]).serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3);
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    HomogeneousIVPSystem dynamics;
    dynamics.vectorfield = eoms;

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0bar(0));
    CHECK(curve.ybar(0, 1) == y0bar(1));
    CHECK(curve.ybar(0, 2) == y0bar(2));
    CHECK(curve.ybar(0, 3) == y0bar(3));

    const size_t nrows = curve.ybar.rows();

    CHECK(std::abs(curve.t(nrows-1) - tf) < TOL_FINE);
    CHECK(std::abs(curve.ybar(nrows-1, 0) - std::exp(m*tf)*y0bar(0)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 1) - std::exp(m*tf)*y0bar(1)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 2) - std::exp(m*tf)*y0bar(2)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 3) - std::exp(m*tf)*y0bar(3)) < 1e-5);
}

TEST_CASE("solve_ivp_1_event", "[integrate]")
{
    /*!
    * Tests solve_ivp against a classical problem with known solution.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;

    const double m = -0.5;
    const double tf = 10.0;
    const double y2cross = 3.5;

    Eigen::VectorXd y0bar(4);
    y0bar(0) = 2.0;
    y0bar(1) = 4.0;
    y0bar(2) = 6.0;
    y0bar(3) = -8.0;

    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    // Do not modify anything below this line.

    const auto eoms = [m](const double t, const CompositeManifold& y)
    {
        const Eigen::VectorXd ybar = std::get<RN>(y.space[0]).serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3);
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    const auto event = [y2cross](const double t, const CompositeManifold& y)
    {
        const Eigen::VectorXd ybar = std::get<RN>(y.space[0]).serialize();
        return ybar(2) - y2cross;
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    HomogeneousIVPSystem dynamics;
    dynamics.vectorfield = eoms;
    dynamics.event = event;

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0bar(0));
    CHECK(curve.ybar(0, 1) == y0bar(1));
    CHECK(curve.ybar(0, 2) == y0bar(2));
    CHECK(curve.ybar(0, 3) == y0bar(3));

    const size_t nrows = curve.ybar.rows();
    const double tcross = std::log(y2cross/y0bar(2))/m;

    CHECK(std::abs(curve.t(nrows-1) - tcross) < TOL_FINE);
    CHECK(std::abs(curve.ybar(nrows-1, 0) - std::exp(m*tcross)*y0bar(0)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 1) - std::exp(m*tcross)*y0bar(1)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 2) - std::exp(m*tcross)*y0bar(2)) < 1e-5);
    CHECK(std::abs(curve.ybar(nrows-1, 3) - std::exp(m*tcross)*y0bar(3)) < 1e-5);
}

TEST_CASE("solve_ivp_1_segmented", "[integrate]")
{
    /*!
    * Tests solve_ivp against a classical problem with known solution.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;
    using namespace Lielab::utils;

    const double m = -0.5;
    const Eigen::VectorXd tspan = linspace(0.0, 10.0, 41);

    Eigen::VectorXd y0bar(4);
    y0bar(0) = 2.0;
    y0bar(1) = 4.0;
    y0bar(2) = 6.0;
    y0bar(3) = -8.0;

    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    // Do not modify anything below this line.

    const auto eoms = [m](const double t, const CompositeManifold& y)
    {
        const Eigen::VectorXd ybar = std::get<RN>(y.space[0]).serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3);
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    HomogeneousIVPSystem dynamics;
    dynamics.vectorfield = eoms;

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    for (size_t ii = 0; ii < tspan.size(); ii++)
    {
        bool val_found = false;
        for (ptrdiff_t jj = 0; jj < curve.t.size(); jj++)
        {
            if (std::abs(tspan(ii) - curve.t(jj)) < 1e-14)
            {
                val_found = true;
            }
        }
        CHECK(val_found);
    }   
}

// Lorenz tests

TEST_CASE("solve_ivp_unwrapped", "[integrate]")
{
    /*!
    * Tests solve_ivp against a function that has not been wrapped.
    */

    using namespace Lielab::integrate;

    const auto eoms = [](const double t, const Eigen::VectorXd& y)
    {
        Eigen::VectorXd dy(3);
        const double sigma = 10.0;
        const double rho = 28.0;
        const double b1 = 8.0;
        const double b2 = 3.0;
        const double beta = b1/b2;

        dy(0) = -beta*y(0) + y(1)*y(2);
        dy(1) = -sigma*y(1) + sigma*y(2);
        dy(2) = -y(0)*y(1) + rho*y(1) - y(2);
        return dy;
    };

    Eigen::VectorXd y0(3);
    y0(0) = 25.0;
    y0(1) = 0.0;
    y0(2) = -20.0;

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = 5.0;

    EuclideanIVPSystem dynamics;
    dynamics.vectorfield = eoms;

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0(0));
    CHECK(curve.ybar(0, 1) == y0(1));
    CHECK(curve.ybar(0, 2) == y0(2));

    const size_t nrows = curve.ybar.rows();

    CHECK(std::abs(curve.t(nrows-1) - 5.0) < TOL_FINE);
    CHECK(std::abs(curve.ybar(nrows-1, 0) - 15.230) < 1e-1);
    CHECK(std::abs(curve.ybar(nrows-1, 1) + 0.797) < 1e-1);
    CHECK(std::abs(curve.ybar(nrows-1, 2) + 1.473) < 1e-1);
}

TEST_CASE("solve_ivp_RNxRN_RN", "[integrate]")
{
    /*!
    * Tests solve_ivp against a function that has been wrapped.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;

    const auto eoms = [](const double t, const CompositeManifold& y)
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

    HomogeneousIVPSystem dynamics;
    dynamics.vectorfield = eoms;

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0bar(0));
    CHECK(curve.ybar(0, 1) == y0bar(1));
    CHECK(curve.ybar(0, 2) == y0bar(2));

    const size_t nrows = curve.ybar.rows();

    CHECK(std::abs(curve.t(nrows-1) - 5.0) < TOL_FINE);
    CHECK(std::abs(curve.ybar(nrows-1, 0) - 15.230) < 1e-1);
    CHECK(std::abs(curve.ybar(nrows-1, 1) + 0.797) < 1e-1);
    CHECK(std::abs(curve.ybar(nrows-1, 2) + 1.473) < 1e-1);
}

TEST_CASE("solve_ivp_GLxRN_RN", "[integrate]")
{
    /*!
    * Tests solve_ivp against a function with custom action GLR x RN -> RN.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;

    const auto eoms = [](const double t, const CompositeManifold& y)
    {
        const Eigen::VectorXd ybar = std::get<RN>(y.space[0]).serialize();

        const double sigma = 10.0;
        const double rho = 28.0;
        const double b1 = 8.0;
        const double b2 = 3.0;
        const double beta = b1/b2;

        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3,3);

        A(0, 0) = -beta;
        A(0, 1) = 0.0;
        A(0, 2) = ybar(1);
        A(1, 0) = 0.0;
        A(1, 1) = -sigma;
        A(1, 2) = sigma;
        A(2, 0) = -ybar(1);
        A(2, 1) = rho;
        A(2, 2) = -1.0;
        
        return CompositeAlgebra({glr(A)});
    };

    const auto action = [](const CompositeGroup& g, const CompositeManifold& y)
    {
        const Eigen::MatrixXd ghat = std::get<GLR>(g.space[0]).get_matrix();
        const Eigen::VectorXd ybar = std::get<RN>(y.space[0]).serialize();
        return CompositeManifold({RN::from_vector(ghat*ybar)});
    };

    Eigen::VectorXd y0bar(3);
    y0bar(0) = 25.0;
    y0bar(1) = 0.0;
    y0bar(2) = -20.0;

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = 5.0;

    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    HomogeneousIVPSystem dynamics;
    dynamics.vectorfield = eoms;
    dynamics.action = action;

    const ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0bar(0));
    CHECK(curve.ybar(0, 1) == y0bar(1));
    CHECK(curve.ybar(0, 2) == y0bar(2));

    const size_t nrows = curve.ybar.rows();

    CHECK(std::abs(curve.t(nrows-1) - 5.0) < TOL_FINE);
    // CHECK(std::abs(curve.ybar(nrows-1, 0) - 15.230) < 1e-1); // TODO:
    // CHECK(std::abs(curve.ybar(nrows-1, 1) + 0.797) < 1e-1);
    // CHECK(std::abs(curve.ybar(nrows-1, 2) + 1.473) < 1e-1);
}

TEST_CASE("solve_ivp_composite1", "[integrate]")
{
    /*!
    * Tests solve_ivp against a function with custom action (SE x GLR) x (SE x RN) -> (SE x RN).
    */

    using namespace Lielab::domain;
    using namespace Lielab::functions;
    using namespace Lielab::integrate;

    const auto eoms = [](const double t, const CompositeManifold& y)
    {
        const double V = 1.0;
        const Eigen::VectorXd lambdabar = std::get<RN>(y.space[1]).serialize();
        const double u = -lambdabar(2);

        const se dx = se::from_vector({V, 0.0, u});
        const glr dlambdabar = -coad(dx);
        
        return CompositeAlgebra({dx, dlambdabar});
    };

    const auto action = [](const CompositeGroup& g, const CompositeManifold& y)
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

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = 1.0;

    const CompositeManifold y0 = CompositeManifold({exp(x0), RN::from_vector(lambda0bar)});

    IVPOptions options;
    // MuntheKaasZanna method(Coefficients::RKDP54_7M); // TODO: Switch to RKV87r (or other default) once better error control is developed
    HomogeneousIVPSystem dynamics;
    dynamics.vectorfield = eoms;
    dynamics.action = action;

    const ODESolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 1) == -1.0);
    CHECK(curve.ybar(0, 3) == 1.0);
    CHECK(curve.ybar(0, 2) == x0bar(0));
    CHECK(curve.ybar(0, 5) == x0bar(1));
    CHECK(curve.ybar(0, 9) == lambda0bar(0));
    CHECK(curve.ybar(0, 10) == lambda0bar(1));
    CHECK(curve.ybar(0, 11) == lambda0bar(2));

    const size_t nrows = curve.ybar.rows();

    CHECK(std::abs(curve.t(nrows-1) - 1.0) < TOL_FINE);
    // CHECK(std::abs(curve.ybar(nrows-1, 2) - 0.12732395447351627) < 1e-1); // TODO:
    // CHECK(std::abs(curve.ybar(nrows-1, 5) + 0.0) < 1e-1);
    // CHECK(std::abs(curve.ybar(nrows-1, 9) - lambda0bar(0)) < 1e-1);
    // CHECK(std::abs(curve.ybar(nrows-1, 10) + lambda0bar(1)) < 1e-1);
    // CHECK(std::abs(curve.ybar(nrows-1, 11) - lambda0bar(2)) < 1e-1);
}

// Technical tests

// Lielab::integrate::vectorfield_t vfex1()
// {
//     /*!
//     * Constructs and returns eom vfex1.
//     * 
//     * Engø, Kenth, Arne Marthinsen, and Hans Z. Munthe-Kaas. "DiffMan: An object-oriented MATLAB
//     * toolbox for solving differential equations on manifolds." Applied numerical mathematics
//     * 39.3-4 (2001): 323-347.
//     */

//     using namespace Lielab::domain;
//     using namespace Lielab::integrate;

//     vectorfield_t out = [](const double t, const CompositeGroup& y)
//     {
//         const SO x = std::get<SO>(y.space[0]);
//         const so dx = so::from_vector({std::pow(t, 2), 1.0, -t});
//         return CompositeAlgebra{dx};
//     };

//     return out;
// }

// std::function<Lielab::domain::CompositeAlgebra(double, Lielab::domain::CompositeGroup)> vfex2()
// {
//     /*!
//      * Constructs and returns eom vfex2.
//      *
//      * Engø, Kenth, Arne Marthinsen, and Hans Z. Munthe-Kaas. "DiffMan: An object-oriented MATLAB
//      * toolbox for solving differential equations on manifolds." Applied numerical mathematics
//      * 39.3-4 (2001): 323-347.
//      */

//     using namespace Lielab::domain;
//     using namespace Lielab::integrate;

//     vectorfield_t out = [](const double t, const CompositeGroup& y)
//     {
//         const RN x = std::get<RN>(y.space[0]);

//         const double sigma = 10.0;
//         const double rho = 28.0;
//         const double b1 = 8.0;
//         const double b2 = 3.0;
//         const double beta = b1/b2;

//         rn dx(4);
//         dx(0) = -beta*x(0) + x(1)*x(2);
//         dx(1) = -sigma*x(1) + sigma*x(2);
//         dx(2) = -x(0)*x(1) + rho*x(1) - x(2);

//         return CompositeAlgebra{dx};
//     };

//     return out;
// }


// TEST_CASE("MuntheKaas_vfex1", "[topos]")
// {
//     /*!
//     * Tests the MuntheKaas function with vfex1.
//     */

//     Lielab::domain::SO y0(3);
//     Lielab::integrate::MuntheKaas ts;
//     std::function<Lielab::domain::CompositeAlgebra(double, Lielab::domain::CompositeGroup)> vf = vfex1();

//     Lielab::domain::CompositeGroup M0{y0};

//     Lielab::integrate::TSOutput out = ts(vf, M0, 1.0, 0.02);

//     Lielab::domain::SO y1 = std::get<Lielab::domain::SO>(out.low.space[0]);

//     CHECK(std::abs(y1(0,0) - 0.999596034819844) < TOL_FINE);
//     CHECK(std::abs(y1(0,1) - 0.020398551422611) < TOL_FINE);
//     CHECK(std::abs(y1(0,2) - 0.019790560181659) < TOL_FINE);
//     CHECK(std::abs(y1(1,0) + 0.019990512514326) < TOL_FINE);
//     CHECK(std::abs(y1(1,1) - 0.999587901244249) < TOL_FINE);
//     CHECK(std::abs(y1(1,2) + 0.020601143063711) < TOL_FINE);
//     CHECK(std::abs(y1(2,0) + 0.020202637992582) < TOL_FINE);
//     CHECK(std::abs(y1(2,1) - 0.020197197478265) < TOL_FINE);
//     CHECK(std::abs(y1(2,2) - 0.999591880035129) < TOL_FINE);
// }

// TEST_CASE("MuntheKaas_vfex2", "[topos]")
// {
//     /*!
//     * Tests the MuntheKaas function with vfex2.
//     */

//     Lielab::domain::RN y0(4);
//     Lielab::integrate::MuntheKaas ts;
//     std::function<Lielab::domain::CompositeAlgebra(double, Lielab::domain::CompositeGroup)> vf = vfex2();

//     y0(0) = 25.0;
//     y0(1) = 0.0;
//     y0(2) = -20.0;

//     Lielab::domain::CompositeGroup M0{y0};

//     Lielab::integrate::TSOutput out = ts(vf, M0, 0.0, 0.02);

//     Lielab::domain::RN y1 = std::get<Lielab::domain::RN>(out.low.space[0]);
    
//     CHECK(std::abs(y1(0) - 24.425986197956878) < TOL_FINE);
//     CHECK(std::abs(y1(1) + 3.596428324678375) < TOL_FINE);
//     CHECK(std::abs(y1(2) + 19.733395791914329) < TOL_FINE);
// }

// TEST_CASE("Flow_fails", "[topos]")
// {
//     /*!
//     * Test cases where Flow should fail.
//     */

//     Lielab::domain::SO y0(3);
//     Lielab::domain::CompositeGroup M0{y0};
//     std::function<Lielab::domain::CompositeAlgebra(double, Lielab::domain::CompositeGroup)> vf = vfex1();
//     std::vector<double> tspan;
//     Lielab::integrate::Flow f;

//     // tspan has no values
//     CHECK_THROWS(f(vf, tspan, M0));

//     // tspan has 1 value
//     tspan.push_back(0.0);
//     CHECK_THROWS(f(vf, tspan, M0));
// }

// TEST_CASE("solve_ivp_copy_output", "[integrate]")
// {
//     /*!
//     * Tests that Flows outputs are copied.
//     */

//     using namespace Lielab::integrate;

//     Flow f;

//     std::function<Eigen::VectorXd(const double, const Eigen::VectorXd&)> vf = [](const double t, const Eigen::VectorXd& y)
//     {
//         Eigen::VectorXd out(2);
//         out(0) = y(1);
//         out(1) = -y(0);
//         return out;
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(5.5);

//     Eigen::VectorXd y0_1(2);
//     Eigen::VectorXd y0_2(2);

//     y0_1 << 1.0, 0.0;
//     y0_2 << 2.0, 0.0;

//     ODESolution curve1 = solve_ivp(vf, tspan, y0_1);
//     ODESolution curve2 = solve_ivp(vf, tspan, y0_2);

//     const size_t L1 = curve1.t.size();
//     const size_t L2 = curve2.t.size();

//     CHECK(std::abs(curve1.ybar(L1-1, 0) - curve2.ybar(L2-1, 0)) >= 1e-4);
//     CHECK(std::abs(curve1.ybar(L1-1, 1) - curve2.ybar(L2-1, 1)) >= 1e-4);

// }

// TEST_CASE("Flow_vfex2_rk45_fixed", "[topos]")
// {
//     Lielab::domain::RN y0(4);
//     std::function<Lielab::domain::CompositeAlgebra(double, Lielab::domain::CompositeGroup)> vf = vfex2();

//     y0(0) = 25.0;
//     y0(1) = 0.0;
//     y0(2) = -20.0;

//     Lielab::domain::CompositeGroup M0{y0};

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(5.0);

//     Lielab::integrate::Flow f;
//     f.variable_time_step = false;
//     f.dt = 0.02;

//     Lielab::integrate::IntegralCurve curve = f(vf, tspan, M0);

//     CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
//     CHECK(curve.y(0, 0) == y0(0));
//     CHECK(curve.y(0, 1) == y0(1));
//     CHECK(curve.y(0, 2) == y0(2));

//     size_t nrows = curve.y.rows();

//     CHECK(std::abs(curve.t(nrows-1) - 5.0) < TOL_FINE);
//     CHECK(std::abs(curve.y(nrows-1, 0) - 15.210570567999987) < TOL_FINE);
//     CHECK(std::abs(curve.y(nrows-1, 1) + 0.788689660918195) < TOL_FINE);
//     CHECK(std::abs(curve.y(nrows-1, 2) + 1.459476938449221) < TOL_FINE);
// }

// TEST_CASE("Flow_vfex2_rk45_variable", "[topos]")
// {
//     Lielab::domain::RN y0(4);
//     std::function<Lielab::domain::CompositeAlgebra(double, Lielab::domain::CompositeGroup)> vf = vfex2();

//     y0(0) = 25.0;
//     y0(1) = 0.0;
//     y0(2) = -20.0;

//     Lielab::domain::CompositeGroup M0{y0};

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(5.0);

//     Lielab::integrate::Flow f;
//     f.variable_time_step = true;
//     f.dt = 0.02;
//     Lielab::integrate::IntegralCurve curve = f(vf, tspan, M0);

//     CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
//     CHECK(std::abs(curve.t(1) - 0.007703769593747) < TOL_COARSE);
//     CHECK(std::abs(curve.t(2) - 0.015420629134474) < TOL_COARSE);
//     CHECK(std::abs(curve.t(3) - 0.023255332563845) < TOL_COARSE);
//     CHECK(std::abs(curve.t(4) - 0.031246577586516) < TOL_COARSE);

//     CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
//     CHECK(curve.y(0, 0) == y0(0));
//     CHECK(curve.y(0, 1) == y0(1));
//     CHECK(curve.y(0, 2) == y0(2));

//     size_t nrows = curve.y.rows();

//     CHECK(std::abs(curve.t(nrows-1) - 5.0) < TOL_FINE);
//     CHECK(std::abs(curve.y(nrows-1, 0) - 15.230102737555342) < TOL_COARSE);
//     CHECK(std::abs(curve.y(nrows-1, 1) + 0.796697875936802) < TOL_COARSE);
//     CHECK(std::abs(curve.y(nrows-1, 2) + 1.472989006310112) < TOL_COARSE);
// }

// TEST_CASE("A1_vector", "[integrate]")
// {
//     /*!
//     * Problem A-1
//     *
//     * Uses VectorXd representation and events.
//     * 
//     * Source: Thompson, S. A collection of test problems for ordinary differential
//     * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
//     * Oak Ridge National Lab., TN (USA), 1987.
//     */

//     using namespace Lielab::integrate;

//     IVPOptions options;

//     std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vf = [](const double t, const Eigen::VectorXd& y)
//     {
//         const double m = 1.0/4.0;
//         const double w = 8.0;
//         const double k = 2.0;
//         Eigen::VectorXd dy = Eigen::VectorXd::Zero(2);
//         dy << y(1), (w - k*y(1))/m;
//         return dy;
//     };

//     std::function<double(double, Eigen::VectorXd)> vf_event = [](const double t, const Eigen::VectorXd& y)
//     {
//         const double H = 10.0;
//         return H - y(0);
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(5.0);

//     Eigen::VectorXd y0 = Eigen::VectorXd::Zero(2);
//     y0 << 0.0, 0.0;

//     ODESolution out = solve_ivp(vf, tspan, y0, options, vf_event);

//     const size_t L = out.t.size();

//     CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
//     CHECK(std::abs(out.t(L-1) - 2.62499999990522) <= 1.0e-9);

//     for (size_t ii = 0; ii < L; ii++)
//     {
//         CHECK(std::abs(out.ybar(ii, 0) - (4.0*(out.t(ii) + 1.0/8.0*std::exp(-8.0*out.t(ii)) - 1.0/8.0))) <= 1e-6);
//         CHECK(std::abs(out.ybar(ii, 1) - (4.0*(1.0 - std::exp(-8.0*out.t(ii))))) <= 1e-6);
//     }
// }

// TEST_CASE("A1_hom", "[integrate]")
// {
//     /*!
//     * Problem A-1
//     *
//     * Uses hom representation and events.
//     * 
//     * Source: Thompson, S. A collection of test problems for ordinary differential
//     * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
//     * Oak Ridge National Lab., TN (USA), 1987.
//     */

//     using namespace Lielab::domain;
//     using namespace Lielab::integrate;

//     IVPOptions options;

//     vectorfield_t vf = [](const double t, const CompositeGroup& y)
//     {
//         const double m = 1.0/4.0;
//         const double w = 8.0;
//         const double k = 2.0;
//         const RN x = std::get<RN>(y.space[0]);

//         const rn dx = rn::from_vector({x(1), (w - k*x(1))/m});
//         return CompositeAlgebra{dx};
//     };

//     event_t vf_event = [](const double t, const CompositeGroup& y)
//     {
//         const double H = 10.0;
//         const RN x = std::get<RN>(y.space[0]);
//         return H - x(0);
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(5.0);

//     const RN a1y0{0.0, 0.0};
//     const CompositeGroup a1M0{a1y0};

//     ODESolution out = solve_ivp(vf, tspan, a1M0, options, vf_event);

//     const size_t L = out.t.size();

//     CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
//     CHECK(std::abs(out.t(L-1) - 2.62499999990522) <= 1.0e-9);

//     for (size_t ii = 0; ii < L; ii++)
//     {
//         CHECK(std::abs(out.ybar(ii, 0) - (4.0*(out.t(ii) + 1.0/8.0*std::exp(-8.0*out.t(ii)) - 1.0/8.0))) <= 1e-6);
//         CHECK(std::abs(out.ybar(ii, 1) - (4.0*(1.0 - std::exp(-8.0*out.t(ii))))) <= 1e-6);
//     }
// }

// TEST_CASE("A2_vector", "[integrate]")
// {
//     /*!
//     * Problem A-2
//     *
//     * Uses VectorXd representation and events.
//     * 
//     * Source: Thompson, S. A collection of test problems for ordinary differential
//     * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
//     * Oak Ridge National Lab., TN (USA), 1987.
//     */



//     Lielab::integrate::Flow f;
//     f.dt_max = 10000.0;

//     const double g = 32.0;
//     const double R = 4000.0 * 5280.0;
//     const double H = 237000.0 * 5280.0;

//     std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vf = [g, R, H](const double t, const Eigen::VectorXd& y)
//     {
//         Eigen::VectorXd dy = Eigen::VectorXd::Zero(1);
//         dy(0) = -std::sqrt(2.0*g*std::pow(R, 2.0)) * std::sqrt((H - y(0))/(H*y(0)));
//         return dy;
//     };

//     std::function<double(double, Eigen::VectorXd)> vf_event = [g, R, H](const double t, const Eigen::VectorXd& y)
//     {
//         return y(0) - R;
//     };

//     std::function<double(double)> h_to_t = [g, R, H](const double h)
//     {
//         return (std::pow(H, 3.0/2.0) / (8.0*h))*(std::sqrt(h/H - std::pow(h/H, 2.0)) + 1.0/2.0*std::acos(2.0*h/H - 1.0));
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(800000.0);

//     Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
//     y0 << H - 1.0e-6;
//     Lielab::integrate::ODESolution out = f(vf, tspan, y0, vf_event);

//     const size_t L = out.t.size();

//     CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
//     CHECK(std::abs(out.t(L-1) - h_to_t(out.ybar(L-1, 0))) <= 1.0e-1);
// }

// TEST_CASE("A2_hom", "[integrate]")
// {
//     /*!
//     * Problem A-2
//     *
//     * Uses hom representation and events.
//     * 
//     * Source: Thompson, S. A collection of test problems for ordinary differential
//     * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
//     * Oak Ridge National Lab., TN (USA), 1987.
//     */

//     using namespace Lielab::domain;
//     using namespace Lielab::integrate;

//     IVPOptions options;
//     options.dt_max = 10000.0;

//     const double g = 32.0;
//     const double R = 4000.0 * 5280.0;
//     const double H = 237000.0 * 5280.0;

//     vectorfield_t vf = [g, R, H](const double t, const CompositeGroup& y)
//     {
//         const RN x = std::get<RN>(y.space[0]);
//         const rn dx = rn::from_vector({-std::sqrt(2.0*g*std::pow(R, 2.0)) * std::sqrt((H - x(0))/(H*x(0)))});
//         return CompositeAlgebra{dx};
//     };

//     event_t vf_event = [g, R, H](const double t, const CompositeGroup& y)
//     {
//         const RN x = std::get<RN>(y.space[0]);
//         return x(0) - R;
//     };

//     std::function<double(double)> h_to_t = [g, R, H](const double h)
//     {
//         return (std::pow(H, 3.0/2.0) / (8.0*h))*(std::sqrt(h/H - std::pow(h/H, 2.0)) + 1.0/2.0*std::acos(2.0*h/H - 1.0));
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(800000.0);

//     const RN y0 = RN::from_vector({H - 1.0e-6});
//     const ODESolution out = solve_ivp(vf, tspan, CompositeGroup{y0}, options, vf_event);

//     const size_t L = out.t.size();

//     CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
//     CHECK(std::abs(out.t(L-1) - h_to_t(out.ybar(L-1, 0))) <= 1.0e-1);
// }

// TEST_CASE("A4_vector", "[integrate]")
// {
//     /*!
//     * Problem A-4
//     *
//     * Uses VectorXd representation and events.
//     * 
//     * Source: Thompson, S. A collection of test problems for ordinary differential
//     * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
//     * Oak Ridge National Lab., TN (USA), 1987.
//     */

//     using namespace Lielab::domain;
//     using namespace Lielab::integrate;

//     constexpr double PI = std::numbers::pi_v<double>;

//     std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> vf = [&PI](const double t, const Eigen::VectorXd& y)
//     {
//         Eigen::VectorXd dy = Eigen::VectorXd::Zero(2);
//         dy << y(1), -(16.0 * std::pow(PI, 2) * std::exp(-2.0*t) - 1.0/4.0)*y(0);
//         return dy;
//     };

//     std::function<double(double)> t_to_y1 = [&PI](const double t)
//     {
//         return std::exp(t/2.0)*std::cos(4.0*PI*std::exp(-t));
//     };

//     std::function<double(double)> t_to_y2 = [&PI](const double t)
//     {
//         return std::exp(t/2.0)*(4.0*PI*std::exp(-t)*std::sin(4.0*PI*std::exp(-t)) + 1.0/2.0*std::cos(4.0*PI*std::exp(-t)));
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(std::log(8.0) - std::log(1.0)); // Root at k=1

//     Eigen::VectorXd y0 = Eigen::VectorXd::Zero(2);
//     y0 << 1.0, 0.5;
    
//     const ODESolution out = solve_ivp(vf, tspan, y0);

//     const size_t L = out.t.size();

//     CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
//     CHECK(std::abs(out.ybar(L-1, 0) - 0.0) <= 1.0e-5);

//     for (size_t ii = 0; ii < L; ii++)
//     {
//         CHECK(std::abs(out.ybar(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-5);
//         CHECK(std::abs(out.ybar(ii, 1) - t_to_y2(out.t(ii))) <= 1.0e-5);
//     }
// }


// TEST_CASE("A4_hom", "[integrate]")
// {
//     /*!
//     * Problem A-4
//     *
//     * Uses hom representation and events.
//     * 
//     * Source: Thompson, S. A collection of test problems for ordinary differential
//     * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
//     * Oak Ridge National Lab., TN (USA), 1987.
//     */

//     using namespace Lielab::domain;
//     using namespace Lielab::integrate;

//     constexpr double PI = std::numbers::pi_v<double>;

//     vectorfield_t vf = [&PI](const double t, const CompositeGroup& y)
//     {
//         const RN x = std::get<RN>(y.space[0]);
//         const rn dx = rn::from_vector({x(1), -(16.0 * std::pow(PI, 2) * std::exp(-2.0*t) - 1.0/4.0)*x(0)});
//         return Lielab::domain::CompositeAlgebra{dx};
//     };

//     std::function<double(double)> t_to_y1 = [&PI](const double t)
//     {
//         return std::exp(t/2.0)*std::cos(4.0*PI*std::exp(-t));
//     };

//     std::function<double(double)> t_to_y2 = [&PI](const double t)
//     {
//         return std::exp(t/2.0)*(4.0*PI*std::exp(-t)*std::sin(4.0*PI*std::exp(-t)) + 1.0/2.0*std::cos(4.0*PI*std::exp(-t)));
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(std::log(8.0) - std::log(1.0)); // Root at k=1

//     const RN y0 = RN::from_vector({1.0, 0.5});
    
//     Lielab::integrate::ODESolution out = solve_ivp(vf, tspan, CompositeGroup{y0});

//     const size_t L = out.t.size();

//     CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
//     CHECK(std::abs(out.ybar(L-1, 0) - 0.0) <= 1.0e-5);

//     for (size_t ii = 0; ii < L; ii++)
//     {
//         CHECK(std::abs(out.ybar(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-5);
//         CHECK(std::abs(out.ybar(ii, 1) - t_to_y2(out.t(ii))) <= 1.0e-5);
//     }
// }

// TEST_CASE("A5_vector", "[integrate]")
// {
//     /*!
//     * Problem A-5
//     *
//     * Uses VectorXd representation and events.
//     * 
//     * Source: Thompson, S. A collection of test problems for ordinary differential
//     * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
//     * Oak Ridge National Lab., TN (USA), 1987.
//     */

//     using namespace Lielab::domain;
//     using namespace Lielab::integrate;

//     IVPOptions options;

//     constexpr double PI = std::numbers::pi_v<double>;

//     std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vf = [&PI](const double t, const Eigen::VectorXd& y)
//     {
//         Eigen::VectorXd dy = Eigen::VectorXd::Zero(2);
//         dy << y(1), -16.0 * std::cos(PI*t/2.0)*y(1) - (64.0 * std::pow(PI, 2) + 64.0*std::pow(std::cos(PI*t/2.0), 2) - 4.0*PI*std::sin(PI*t/2.0))*y(0);
//         return dy;
//     };

//     std::function<double(double, Eigen::VectorXd)> vf_event = [](const double t, const Eigen::VectorXd & y)
//     {
//         if ((t > 1.07) && (t < 1.30))
//         {
//             return -y(0);
//         }
//         return 1.0;
//     };

//     std::function<double(double)> t_to_y1 = [&PI](const double t)
//     {
//         return std::exp(-16/PI*std::sin(PI*t/2))*std::cos(8*PI*t);
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(2.0);

//     Eigen::VectorXd y0 = Eigen::VectorXd::Zero(2);
//     y0 << 1.0, -8.0;
    
//     const ODESolution out = solve_ivp(vf, tspan, y0, options, vf_event);

//     const size_t L = out.t.size();

//     CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
//     CHECK(std::abs(out.t(L-1) - (2.0*10.0-1.0)/16.0) <= 1.0e-5);

//     for (size_t ii = 0; ii < L; ii++)
//     {
//         CHECK(std::abs(out.ybar(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-5);
//     }
// }

// TEST_CASE("A5_hom", "[integrate]")
// {
//     /*!
//     * Problem A-5
//     *
//     * Uses hom representation and events.
//     * 
//     * Source: Thompson, S. A collection of test problems for ordinary differential
//     * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
//     * Oak Ridge National Lab., TN (USA), 1987.
//     */

//     using namespace Lielab::domain;
//     using namespace Lielab::integrate;

//     IVPOptions options;

//     constexpr double PI = std::numbers::pi_v<double>;

//     vectorfield_t vf = [&PI](const double t, const CompositeGroup& y)
//     {
//         const RN x = std::get<RN>(y.space[0]);
//         const rn dx = rn::from_vector({x(1), -16.0 * std::cos(PI*t/2.0)*x(1) - (64.0 * std::pow(PI, 2) + 64.0*std::pow(std::cos(PI*t/2.0), 2) - 4.0*PI*std::sin(PI*t/2.0))*x(0)});
//         return CompositeAlgebra{dx};
//     };

//     event_t vf_event = [](const double t, const CompositeGroup& y)
//     {
//         const RN x = std::get<RN>(y.space[0]);
//         if ((t > 1.07) && (t < 1.30))
//         {
//             return -x(0);
//         }
//         return 1.0;
//     };

//     std::function<double(double)> t_to_y1 = [&PI](const double t)
//     {
//         return std::exp(-16/PI*std::sin(PI*t/2))*std::cos(8*PI*t);
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(2.0);

//     const RN y0 = RN::from_vector({1.0, -8.0});
    
//     const ODESolution out = solve_ivp(vf, tspan, CompositeGroup{y0}, options, vf_event);

//     const size_t L = out.t.size();

//     CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
//     CHECK(std::abs(out.t(L-1) - (2.0*10.0-1.0)/16.0) <= 1.0e-5);

//     for (size_t ii = 0; ii < L; ii++)
//     {
//         CHECK(std::abs(out.ybar(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-5);
//     }
// }

// TEST_CASE("B1_vector", "[integrate]")
// {
//     /*!
//     * Problem B-1
//     *
//     * Uses VectorXd representation and events.
//     * 
//     * Source: Thompson, S. A collection of test problems for ordinary differential
//     * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
//     * Oak Ridge National Lab., TN (USA), 1987.
//     */

//     using namespace Lielab::integrate;

//     const double K = 5.0/3.0;
//     const double L = 1.0/4.0;

//     std::function<Eigen::VectorXd(const double, const Eigen::VectorXd&)> vf = [K, L](const double t, const Eigen::VectorXd& y)
//     {
//         Eigen::VectorXd dy = Eigen::VectorXd::Zero(2);
//         const double dy0 = y(0);
//         double dy1 = 0.0;

//         const double E = y(0) - y(1);

//         if (E < -L/K)
//         {
//             dy1 = -L;
//         }
//         else if ((-L/K <= E) && (E <= L/K))
//         {
//             dy1 = K*E;
//         }
//         else if (L/K < E)
//         {
//             dy1 = L;
//         }

//         dy << dy0, dy1;
//         return dy;
//     };

//     const double t1 = 0.1569;
//     const double y2t1 = 1.0199;

//     std::function<double(double)> t_to_y1 = [](const double t)
//     {
//         return std::exp(t);
//     };

//     std::function<double(double)> t_to_y2 = [K, L, t1, y2t1](const double t)
//     {
//         if (t < t1)
//         {
//             return K/(K+1)*std::exp(t) + 1/(K+1)*std::exp(-K*t);
//         }
            
//         return L*(t - t1) + y2t1;
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(0.5);

//     Eigen::VectorXd y0 = Eigen::VectorXd::Zero(2);
//     y0 << 1.0, 1.0;
    
//     const ODESolution out = solve_ivp(vf, tspan, y0);

//     const size_t len = out.t.size();

//     CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-14);
//     CHECK(std::abs(out.t(len-1) - 0.5) <= 1.0e-14);

//     for (size_t ii = 0; ii < len; ii++)
//     {
//         CHECK(std::abs(out.ybar(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-7);
//         CHECK(std::abs(out.ybar(ii, 1) - t_to_y2(out.t(ii))) <= 1.0e-4); // Answer given in the document is only good to 1e-4 (see t1 and y2(t1))
//     }
// }

// TEST_CASE("B1_hom", "[integrate]")
// {
//     /*!
//     * Problem B-1
//     *
//     * Uses hom representation and events.
//     * 
//     * Source: Thompson, S. A collection of test problems for ordinary differential
//     * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
//     * Oak Ridge National Lab., TN (USA), 1987.
//     */

//     using namespace Lielab::domain;
//     using namespace Lielab::integrate;

//     std::function<CompositeGroup(const CompositeGroup&, const CompositeGroup&)> action = [](const CompositeGroup& g, const CompositeGroup& y)
//     {
//         const GLR _G0 = std::get<GLR>(g.space[0]);
//         const RN _G1 = std::get<RN>(g.space[1]);
//         const RN _Y0 = std::get<RN>(y.space[0]);
//         const RN _Y1 = std::get<RN>(y.space[1]);
        
//         const RN _Y0next = RN::from_vector(_G0.data*_Y0.data);
//         const RN _Y1next = _G1*_Y1;
//         return CompositeGroup{_Y0next, _Y1next};
//     };

//     Flow f;
//     MuntheKaas method;
//     method.Lie_group_action = action;
//     f.method = std::make_shared<MuntheKaas>(method);

//     const double K = 5.0/3.0;
//     const double L = 1.0/4.0;

//     vectorfield_t vf = [K, L](const double t, const CompositeGroup& y)
//     {
//         RN _y0 = std::get<RN>(y.space[0]);
//         RN _y1 = std::get<RN>(y.space[1]);

//         const double y0 = _y0(0);
//         const double y1 = _y1(0);

//         double dy1 = 0.0;

//         const double E = y0 - y1;

//         if (E < -L/K)
//         {
//             dy1 = -L;
//         }
//         else if ((-L/K <= E) && (E <= L/K))
//         {
//             dy1 = K*E;
//         }
//         else if (L/K < E)
//         {
//             dy1 = L;
//         }

//         return CompositeAlgebra({glr::basis(0,1), dy1*rn::basis(0,2)});
//     };

//     const double t1 = 0.1569;
//     const double y2t1 = 1.0199;

//     std::function<double(double)> t_to_y1 = [](const double t)
//     {
//         return std::exp(t);
//     };

//     std::function<double(double)> t_to_y2 = [K, L, t1, y2t1](const double t)
//     {
//         if (t < t1)
//         {
//             return K/(K+1)*std::exp(t) + 1/(K+1)*std::exp(-K*t);
//         }
            
//         return L*(t - t1) + y2t1;
//     };

//     std::vector<double> tspan;
//     tspan.push_back(0.0);
//     tspan.push_back(0.5);

//     const CompositeGroup y0({RN::from_vector({1.0}), RN::from_vector({1.0})});
    
//     const ODESolution out = f(vf, tspan, y0);

//     const size_t len = out.t.size();

//     CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-14);
//     CHECK(std::abs(out.t(len-1) - 0.5) <= 1.0e-14);

//     for (size_t ii = 0; ii < len; ii++)
//     {
//         CHECK(std::abs(out.ybar(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-14); // Answer is analytic
//         CHECK(std::abs(out.ybar(ii, 1) - t_to_y2(out.t(ii))) <= 1.0e-4); // Answer given in the document is only good to 1e-4 (see t1 and y2(t1))
//     }
// }
