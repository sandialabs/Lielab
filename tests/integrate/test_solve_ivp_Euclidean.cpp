#include <functional>
#include <memory>
#include <numbers>

#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

/*
Unit-like tests for solve_ivp with Euclidean systems.
*/

/*
Error handling
*/

TEST_CASE("solve_ivp_Euclidean_nan", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when nans appear in the vectorfield.
    */
    
    using namespace Lielab::integrate;

    const double m = -0.5;
    const double tf = 10.0;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    const auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
    {
        Eigen::VectorXd dy(4);

        dy(0) = m*y(0);
        dy(1) = m*y(1);
        dy(2) = m*y(2);
        dy(3) = m*y(3)*std::numeric_limits<double>::quiet_NaN();
        return dy;
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    EuclideanIVPSystem dynamics(eoms);

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(curve.status == -1);
    CHECK(curve.success == false);
}

TEST_CASE("solve_ivp_Euclidean_inf", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when infs appear in the vectorfield.
    */
    
    using namespace Lielab::integrate;

    const double m = -0.5;
    const double tf = 10.0;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    const auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
    {
        Eigen::VectorXd dy(4);

        dy(0) = m*y(0);
        dy(1) = m*y(1);
        dy(2) = m*y(2);
        dy(3) = m*y(3)*std::numeric_limits<double>::infinity();
        return dy;
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    EuclideanIVPSystem dynamics(eoms);

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(curve.status == -2);
    CHECK(curve.success == false);
}

TEST_CASE("solve_ivp_Euclidean_tspan_short", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when tspan is too short.
    */
    
    using namespace Lielab::integrate;

    const double m = -0.5;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    const auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
    {
        Eigen::VectorXd dy(4);

        dy(0) = m*y(0);
        dy(1) = m*y(1);
        dy(2) = m*y(2);
        dy(3) = m*y(3);
        return dy;
    };

    Eigen::VectorXd tspan(1);
    tspan(0) = 0.0;

    EuclideanIVPSystem dynamics(eoms);

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(curve.status == -3);
    CHECK(curve.success == false);
}

TEST_CASE("solve_ivp_Euclidean_tspan_repeat", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when tspan has repeated values.
    */
    
    using namespace Lielab::integrate;

    const double m = -0.5;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    const auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
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
    tspan(1) = 0.0;

    EuclideanIVPSystem dynamics(eoms);

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(curve.status == -20);
    CHECK(curve.success == false);
}

TEST_CASE("solve_ivp_Euclidean_tspan_descending", "[integrate]")
{
    /*
    Tests solve_ivp for error handling when tspan is in descending order.
    */
    
    using namespace Lielab::integrate;

    const double m = -0.5;
    const double tf = 10.0;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    const auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
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
    tspan(1) = -tf;

    EuclideanIVPSystem dynamics(eoms);

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(curve.status == -20);
    CHECK(curve.success == false);
}

/*
Basic features
*/

TEST_CASE("solve_ivp_Euclidean", "[integrate]")
{
    /*!
    Tests solve_ivp against a classical problem with known solution.
    */

    using namespace Lielab::integrate;
    using Eigen::placeholders::last;

    const double m = -0.5;
    const double tf = 10.0;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    const auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
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

    EuclideanIVPSystem dynamics(eoms);

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) <1e-15);
    CHECK(curve.ybar(0, 0) == y0(0));
    CHECK(curve.ybar(0, 1) == y0(1));
    CHECK(curve.ybar(0, 2) == y0(2));
    CHECK(curve.ybar(0, 3) == y0(3));

    CHECK(std::abs(curve.t(last) - tf) < 1e-15);
    CHECK(std::abs(curve.ybar(last, 0) - std::exp(m*tf)*y0(0)) < 1e-5);
    CHECK(std::abs(curve.ybar(last, 1) - std::exp(m*tf)*y0(1)) < 1e-5);
    CHECK(std::abs(curve.ybar(last, 2) - std::exp(m*tf)*y0(2)) < 1e-5);
    CHECK(std::abs(curve.ybar(last, 3) - std::exp(m*tf)*y0(3)) < 1e-5);
}

TEST_CASE("solve_ivp_Euclidean_tol", "[integrate]")
{
    /*!
    Tests solve_ivp for reporting tolerance issues when the step cannot adapt properly.
    */

    using namespace Lielab::integrate;

    const double m = -0.5;
    const double tf = 10.0;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    const auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
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

    EuclideanIVPSystem dynamics(eoms);

    IVPOptions options;
    options.dt_min = 5.0;

    ODESolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(curve.status == 4);
    CHECK(curve.success == true);

    CHECK(curve.t.size() == 3);
    // CHECK len y
    CHECK(curve.ybar.rows() == 3);
    CHECK(curve.ybar.cols() == 4);
    // CHECK len theta
    CHECK(curve.thetabar.rows() == 3);
    CHECK(curve.thetabar.cols() == 4);
}

TEST_CASE("solve_ivp_Euclidean_event", "[integrate]")
{
    /*!
    Tests solve_ivp against a classical problem with known solution.
    */

    using namespace Lielab::integrate;
    using Eigen::placeholders::last;

    const double m = -0.5;
    const double tf = 10.0;
    const double y2cross = 3.5;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    // Do not modify anything below this line.

    const auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
    {
        Eigen::VectorXd dy(4);

        dy(0) = m*y(0);
        dy(1) = m*y(1);
        dy(2) = m*y(2);
        dy(3) = m*y(3);
        return dy;
    };

    const auto event = [y2cross](const double t, const Eigen::VectorXd& y) -> double
    {
        return y(2) - y2cross;
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    EuclideanIVPSystem dynamics(eoms);
    dynamics.event = event;

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < 1e-15);
    CHECK(curve.ybar(0, 0) == y0(0));
    CHECK(curve.ybar(0, 1) == y0(1));
    CHECK(curve.ybar(0, 2) == y0(2));
    CHECK(curve.ybar(0, 3) == y0(3));

    const double tcross = std::log(y2cross/y0(2))/m;

    CHECK(std::abs(curve.t(last) - tcross) < 1e-8);
    CHECK(std::abs(curve.ybar(last, 0) - std::exp(m*tcross)*y0(0)) < 1e-5);
    CHECK(std::abs(curve.ybar(last, 1) - std::exp(m*tcross)*y0(1)) < 1e-5);
    CHECK(std::abs(curve.ybar(last, 2) - std::exp(m*tcross)*y0(2)) < 1e-5);
    CHECK(std::abs(curve.ybar(last, 3) - std::exp(m*tcross)*y0(3)) < 1e-5);
}

TEST_CASE("solve_ivp_Euclidean_event_tol", "[integrate]")
{
    /*!
    Tests solve_ivp handling of events and tolerance issues with a classical problem with known solution.
    */

    using namespace Lielab::integrate;

    const double m = -0.5;
    const double tf = 10.0;
    const double y2cross = 0.25;

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    const auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
    {
        Eigen::VectorXd dy(4);

        dy(0) = m*y(0);
        dy(1) = m*y(1);
        dy(2) = m*y(2);
        dy(3) = m*y(3);
        return dy;
    };

    const auto event = [y2cross](const double t, const Eigen::VectorXd& y) -> double
    {
        return y(2) - y2cross;
    };

    Eigen::VectorXd tspan(2);
    tspan(0) = 0.0;
    tspan(1) = tf;

    EuclideanIVPSystem dynamics(eoms);
    dynamics.event = event;

    IVPOptions options;
    options.dt_min = 5.0;

    ODESolution curve = solve_ivp(dynamics, tspan, y0, options);

    CHECK(curve.status == 5);
    CHECK(curve.success == true);

    CHECK(curve.t.size() == 3);
    // CHECK len y
    CHECK(curve.ybar.rows() == 3);
    CHECK(curve.ybar.cols() == 4);
    // CHECK len theta
    CHECK(curve.thetabar.rows() == 3);
    CHECK(curve.thetabar.cols() == 4);
}

TEST_CASE("solve_ivp_Euclidean_segmented", "[integrate]")
{
    /*!
    Tests solve_ivp exactly reporting elements in tspan on a classical problem with known solution.
    */

    using namespace Lielab::integrate;
    using namespace Lielab::utils;

    const double m = -0.5;
    const Eigen::VectorXd tspan = linspace(0.0, 10.0, 41);

    Eigen::VectorXd y0(4);
    y0(0) = 2.0;
    y0(1) = 4.0;
    y0(2) = 6.0;
    y0(3) = -8.0;

    const auto eoms = [m](const double t, const Eigen::VectorXd& y) -> Eigen::VectorXd
    {
        Eigen::VectorXd dy(4);

        dy(0) = m*y(0);
        dy(1) = m*y(1);
        dy(2) = m*y(2);
        dy(3) = m*y(3);
        return dy;
    };

    EuclideanIVPSystem dynamics(eoms);

    ODESolution curve = solve_ivp(dynamics, tspan, y0);

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
