#include <functional>
#include <memory>
#include <numbers>

#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

/*
Unit-like tests for solve_ivp with Homogeneous systems.
*/

/*
Error handling
*/

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_nan", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when nans appear in the vectorfield.
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

    const auto eoms = [m](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3)*std::numeric_limits<double>::quiet_NaN();
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    HomogeneousIVPSystem dynamics(eoms);

    IVPOptions options;
    options.method = IVPMethod::MuntheKaas;

    IVPSolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(curve.status == IVPStatus::ERROR_NANS_IN_VF);
    CHECK(curve.success == false);
}

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_inf", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when infs appear in the vectorfield.
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

    const auto eoms = [m](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3)*std::numeric_limits<double>::infinity();
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    HomogeneousIVPSystem dynamics(eoms);

    IVPOptions options;
    options.method = IVPMethod::MuntheKaas;

    IVPSolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(curve.status == IVPStatus::ERROR_INFS_IN_VF);
    CHECK(curve.success == false);
}

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_tspan_short", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when tspan is too short.
    */
    
    using namespace Lielab::domain;
    using namespace Lielab::integrate;

    const double m = -0.5;

    Eigen::VectorXd y0bar(4);
    y0bar(0) = 2.0;
    y0bar(1) = 4.0;
    y0bar(2) = 6.0;
    y0bar(3) = -8.0;

    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    const auto eoms = [m](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3);
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    Eigen::VectorXd tspan(1);
    tspan(0) = 0.0;

    HomogeneousIVPSystem dynamics(eoms);

    IVPOptions options;
    options.method = IVPMethod::MuntheKaas;

    CHECK_THROWS(solve_ivp(dynamics, tspan, y0, options));
}

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_tspan_repeat", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when tspan has repeated values.
    */
    
    using namespace Lielab::domain;
    using namespace Lielab::integrate;

    const double m = -0.5;

    Eigen::VectorXd y0bar(4);
    y0bar(0) = 2.0;
    y0bar(1) = 4.0;
    y0bar(2) = 6.0;
    y0bar(3) = -8.0;

    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    const auto eoms = [m](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3);
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = 0.0;

    HomogeneousIVPSystem dynamics(eoms);

    IVPOptions options;
    options.method = IVPMethod::MuntheKaas;

    CHECK_THROWS(solve_ivp(dynamics, tspan, y0, options));
}

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_tspan_descending", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when tspan is in descending order.
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

    const auto eoms = [m](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3);
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = -tf;

    HomogeneousIVPSystem dynamics(eoms);

    IVPOptions options;
    options.method = IVPMethod::MuntheKaas;

    CHECK_THROWS(solve_ivp(dynamics, tspan, y0, options));
}

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_topos", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when user provided topology comes out incorrect.
    */
    
    using namespace Lielab::domain;
    using namespace Lielab::integrate;
    using namespace Lielab::utils;

    const auto eoms = [](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        return CompositeAlgebra({so::from_vector({1.0, 2.0, 3.0})});
    };

    const auto action = [](const CompositeGroup& g, const CompositeManifold& y)
    {
        const SO& g0 = g[0];
        const SO& y0 = y[0];
        return CompositeManifold({g0*y0, SU(2)});
    };

    const Eigen::VectorXd tspan = to_VectorXd({0.0, 5.0});

    HomogeneousIVPSystem dynamics(eoms);
    dynamics.action = action;

    IVPOptions options;
    options.method = IVPMethod::MuntheKaas;

    CHECK_THROWS(solve_ivp(dynamics, tspan, CompositeManifold({SO(3)}), options));
}

/*
Basic features
*/

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas", "[integrate]")
{
    /*!
    Tests solve_ivp against a classical problem with known solution.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;
    using Eigen::placeholders::last;

    const double m = -0.5;
    const double tf = 10.0;

    Eigen::VectorXd y0bar(4);
    y0bar(0) = 2.0;
    y0bar(1) = 4.0;
    y0bar(2) = 6.0;
    y0bar(3) = -8.0;

    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    const auto eoms = [m](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
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

    HomogeneousIVPSystem dynamics(eoms);
    
    IVPOptions options;
    options.method = IVPMethod::MuntheKaas;

    IVPSolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(curve.status == IVPStatus::SUCCESS);
    CHECK(curve.success == true);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0bar(0));
    CHECK(curve.ybar(0, 1) == y0bar(1));
    CHECK(curve.ybar(0, 2) == y0bar(2));
    CHECK(curve.ybar(0, 3) == y0bar(3));

    CHECK(std::abs(curve.t(last) - tf) < TOL_FINE);
    CHECK(std::abs(curve.ybar(last, 0) - std::exp(m*tf)*y0bar(0)) < 1e-5);
    CHECK(std::abs(curve.ybar(last, 1) - std::exp(m*tf)*y0bar(1)) < 1e-5);
    CHECK(std::abs(curve.ybar(last, 2) - std::exp(m*tf)*y0bar(2)) < 1e-5);
    CHECK(std::abs(curve.ybar(last, 3) - std::exp(m*tf)*y0bar(3)) < 1e-5);
}

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_tol", "[integrate]")
{
    /*!
    Tests solve_ivp for reporting tolerance issues when the step cannot adapt properly.
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

    const auto eoms = [m](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
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

    HomogeneousIVPSystem dynamics(eoms);

    IVPOptions options;
    options.dt_min = 5.0;
    options.method = IVPMethod::MuntheKaas;

    IVPSolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(curve.status == IVPStatus::SUCCESS_BUT_TOL);
    CHECK(curve.success == true);

    CHECK(curve.t.size() == 3);
    // CHECK len y
    CHECK(curve.ybar.rows() == 3);
    CHECK(curve.ybar.cols() == 4);
    // CHECK len theta
    CHECK(curve.thetabar.rows() == 3);
    CHECK(curve.thetabar.cols() == 4);
}

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_event", "[integrate]")
{
    /*!
    Tests solve_ivp against a classical problem with known solution.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;
    using Eigen::placeholders::last;

    const double m = -0.5;
    const double tf = 10.0;
    const double y2cross = 3.5;

    Eigen::VectorXd y0bar(4);
    y0bar(0) = 2.0;
    y0bar(1) = 4.0;
    y0bar(2) = 6.0;
    y0bar(3) = -8.0;

    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    const auto eoms = [m](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3);
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    const auto event = [y2cross](const double t, const CompositeManifold& y) -> double
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        return ybar(2) - y2cross;
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    HomogeneousIVPSystem dynamics(eoms);
    dynamics.event = event;

    IVPOptions options;
    options.method = IVPMethod::MuntheKaas;

    IVPSolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(curve.status == IVPStatus::SUCCESS_EVENT);
    CHECK(curve.success == true);

    CHECK(std::abs(curve.t(0) - 0.0) < options.reltol*10.0);
    CHECK(curve.ybar(0, 0) == y0bar(0));
    CHECK(curve.ybar(0, 1) == y0bar(1));
    CHECK(curve.ybar(0, 2) == y0bar(2));
    CHECK(curve.ybar(0, 3) == y0bar(3));

    const double tcross = std::log(y2cross/y0bar(2))/m;

    CHECK(std::abs(curve.t(last) - tcross) < options.reltol*10.0);
    CHECK(std::abs(curve.ybar(last, 0) - std::exp(m*tcross)*y0bar(0)) < options.reltol*10.0);
    CHECK(std::abs(curve.ybar(last, 1) - std::exp(m*tcross)*y0bar(1)) < options.reltol*10.0);
    CHECK(std::abs(curve.ybar(last, 2) - std::exp(m*tcross)*y0bar(2)) < options.reltol*10.0);
    CHECK(std::abs(curve.ybar(last, 3) - std::exp(m*tcross)*y0bar(3)) < options.reltol*10.0);
}

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_event_tol", "[integrate]")
{
    /*!
    Tests solve_ivp handling of events and tolerance issues with a classical problem with known solution.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;

    const double m = -0.5;
    const double tf = 10.0;
    const double y2cross = 0.25;

    Eigen::VectorXd y0bar(4);
    y0bar(0) = 2.0;
    y0bar(1) = 4.0;
    y0bar(2) = 6.0;
    y0bar(3) = -8.0;

    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    const auto eoms = [m](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3);
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    const auto event = [y2cross](const double t, const CompositeManifold& y) -> double
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        return ybar(2) - y2cross;
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    HomogeneousIVPSystem dynamics(eoms);
    dynamics.event = event;

    IVPOptions options;
    options.dt_min = 5.0;
    options.method = IVPMethod::MuntheKaas;

    IVPSolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(curve.status == IVPStatus::SUCCESS_EVENT_BUT_TOL);
    CHECK(curve.success == true);

    CHECK(curve.t.size() == 3);
    // CHECK len y
    CHECK(curve.ybar.rows() == 3);
    CHECK(curve.ybar.cols() == 4);
    // CHECK len theta
    CHECK(curve.thetabar.rows() == 3);
    CHECK(curve.thetabar.cols() == 4);
}

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_segmented", "[integrate]")
{
    /*!
    Tests solve_ivp exactly reporting elements in tspan on a classical problem with known solution.
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

    const auto eoms = [m](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();
        Eigen::VectorXd dy(4);

        dy(0) = m*ybar(0);
        dy(1) = m*ybar(1);
        dy(2) = m*ybar(2);
        dy(3) = m*ybar(3);
        return CompositeAlgebra({rn::from_vector(dy)});
    };

    HomogeneousIVPSystem dynamics(eoms);

    IVPOptions options;
    options.method = IVPMethod::MuntheKaas;

    IVPSolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(curve.status == IVPStatus::SUCCESS);
    CHECK(curve.success == true);

    for (ptrdiff_t ii = 0; ii < tspan.size(); ii++)
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

/*
Advanced features
*/

TEST_CASE("solve_ivp_Homogeneous_MuntheKaas_Lorenz_GLRxRN_RN", "[integrate]")
{
    /*!
    * Tests solve_ivp against a function with custom action GLR x RN -> RN.
    */

    using namespace Lielab::domain;
    using namespace Lielab::integrate;
    using namespace Lielab::utils;
    using Eigen::placeholders::last;

    const auto eoms = [](const double t, const CompositeManifold& y) -> CompositeAlgebra
    {
        const RN& y0 = y[0];
        const Eigen::VectorXd ybar = y0.serialize();

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

    const auto action = [](const CompositeGroup& g, const CompositeManifold& y) -> CompositeManifold
    {
        const GLR& g0 = g[0];
        const RN& y0 = y[0];
        const Eigen::MatrixXd ghat = g0.get_matrix();
        const Eigen::VectorXd ybar = y0.serialize();
        return CompositeManifold({RN::from_vector(ghat*ybar)});
    };

    Eigen::VectorXd y0bar(3);
    y0bar(0) = 25.0;
    y0bar(1) = 0.0;
    y0bar(2) = -20.0;

    Eigen::VectorXd tspan = linspace<double>(0.0, 5.0, 2);

    const CompositeManifold y0 = CompositeManifold({RN::from_vector(y0bar)});

    HomogeneousIVPSystem dynamics(eoms);
    dynamics.action = action;

    IVPOptions options;
    options.method = IVPMethod::MuntheKaas;

    const IVPSolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.ybar(0, 0) == y0bar(0));
    CHECK(curve.ybar(0, 1) == y0bar(1));
    CHECK(curve.ybar(0, 2) == y0bar(2));

    CHECK(std::abs(curve.t(last) - 5.0) < TOL_FINE);
    CHECK(std::abs(curve.ybar(last, 0) - 15.230) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 1) + 0.797) < 1e-1);
    CHECK(std::abs(curve.ybar(last, 2) + 1.473) < 1e-1);
}
