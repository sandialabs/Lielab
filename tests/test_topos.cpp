#include <functional>

// Lorenz EOMs
Eigen::VectorXd eoms_lorenz_unwrapped(double t, Eigen::VectorXd X)
{
    Eigen::VectorXd dx(3);
    double sigma = 10.0;
    double rho = 28.0;
    double b1 = 8.0;
    double b2 = 3.0;
    double beta = b1/b2;

    dx(0) = -beta*X(0) + X(1)*X(2);
    dx(1) = -sigma*X(1) + sigma*X(2);
    dx(2) = -X(0)*X(1) + rho*X(1) - X(2);
    return dx;
}

TEST_CASE("MuntheKaas_vfex1", "[topos]")
{
    /*!
    * Tests the MuntheKaas function with vfex1.
    */

    lielab::domain::SO y0(3);
    lielab::topos::MuntheKaas ts;
    lielab::dynamics::vectorfield vf = lielab::dynamics::vfex1();

    lielab::domain::hmlie M0{y0};

    lielab::topos::TSOutput out = ts(vf, M0, 1.0, 0.02);

    lielab::domain::SO y1 = std::get<lielab::domain::SO>(out.low.space[0]);

    CHECK(std::abs(y1(0,0) - 0.999596034819844) < TOL_FINE);
    CHECK(std::abs(y1(0,1) - 0.020398551422611) < TOL_FINE);
    CHECK(std::abs(y1(0,2) - 0.019790560181659) < TOL_FINE);
    CHECK(std::abs(y1(1,0) + 0.019990512514326) < TOL_FINE);
    CHECK(std::abs(y1(1,1) - 0.999587901244249) < TOL_FINE);
    CHECK(std::abs(y1(1,2) + 0.020601143063711) < TOL_FINE);
    CHECK(std::abs(y1(2,0) + 0.020202637992582) < TOL_FINE);
    CHECK(std::abs(y1(2,1) - 0.020197197478265) < TOL_FINE);
    CHECK(std::abs(y1(2,2) - 0.999591880035129) < TOL_FINE);
}

TEST_CASE("MuntheKaas_vfex2", "[topos]")
{
    /*!
    * Tests the MuntheKaas function with vfex2.
    */

    lielab::domain::RN y0(4);
    lielab::topos::MuntheKaas ts;
    lielab::dynamics::vectorfield vf = lielab::dynamics::vfex2();

    y0(0) = 25.0;
    y0(1) = 0.0;
    y0(2) = -20.0;

    lielab::domain::hmlie M0{y0};

    lielab::topos::TSOutput out = ts(vf, M0, 0.0, 0.02);

    lielab::domain::RN y1 = std::get<lielab::domain::RN>(out.low.space[0]);
    
    CHECK(std::abs(y1(0) - 24.425986197956878) < TOL_FINE);
    CHECK(std::abs(y1(1) + 3.596428324678375) < TOL_FINE);
    CHECK(std::abs(y1(2) + 19.733395791914329) < TOL_FINE);
}

TEST_CASE("Flow_fails", "[topos]")
{
    /*!
    * Test cases where Flow should fail.
    */

    lielab::domain::SO y0(3);
    lielab::domain::hmlie M0{y0};
    lielab::dynamics::vectorfield vf = lielab::dynamics::vfex1();
    std::vector<double> tspan;
    lielab::topos::Flow f;

    // tspan has no values
    CHECK_THROWS(f(vf, tspan, M0));

    // tspan has 1 value
    tspan.push_back(0.0);
    CHECK_THROWS(f(vf, tspan, M0));
}

TEST_CASE("Flow_copy_output", "[topos]")
{
    /*!
    * Tests that Flows outputs are copied.
    */

    lielab::topos::Flow f;

    std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vf = [](const double t, const Eigen::VectorXd & y)
    {
        Eigen::VectorXd out(2);
        out << y(1), -y(0);
        return out;
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(5.5);

    Eigen::VectorXd y0_1(2);
    Eigen::VectorXd y0_2(2);

    y0_1 << 1.0, 0.0;
    y0_2 << 2.0, 0.0;

    lielab::topos::IntegralCurve curve1 = f(vf, tspan, y0_1);
    lielab::topos::IntegralCurve curve2 = f(vf, tspan, y0_2);

    const size_t L1 = curve1.t.size();
    const size_t L2 = curve2.t.size();

    CHECK(std::abs(curve1.y(L1-1, 0) - curve2.y(L2-1, 0)) >= 1e-4);
    CHECK(std::abs(curve1.y(L1-1, 1) - curve2.y(L2-1, 1)) >= 1e-4);

}

TEST_CASE("Flow_vfex2_rk45_fixed", "[topos]")
{
    lielab::domain::RN y0(4);
    lielab::dynamics::vectorfield vf = lielab::dynamics::vfex2();

    y0(0) = 25.0;
    y0(1) = 0.0;
    y0(2) = -20.0;

    lielab::domain::hmlie M0{y0};

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(5.0);

    lielab::topos::Flow f;
    f.variable_time_step = false;
    f.dt = 0.02;

    lielab::topos::IntegralCurve curve = f(vf, tspan, M0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.y(0, 0) == y0(0));
    CHECK(curve.y(0, 1) == y0(1));
    CHECK(curve.y(0, 2) == y0(2));

    size_t nrows = curve.y.rows();

    CHECK(std::abs(curve.t(nrows-1) - 5.0) < TOL_FINE);
    CHECK(std::abs(curve.y(nrows-1, 0) - 15.210570567999987) < TOL_FINE);
    CHECK(std::abs(curve.y(nrows-1, 1) + 0.788689660918195) < TOL_FINE);
    CHECK(std::abs(curve.y(nrows-1, 2) + 1.459476938449221) < TOL_FINE);
}

TEST_CASE("Flow_vfex2_rk45_variable", "[topos]")
{
    lielab::domain::RN y0(4);
    lielab::dynamics::vectorfield vf = lielab::dynamics::vfex2();

    y0(0) = 25.0;
    y0(1) = 0.0;
    y0(2) = -20.0;

    lielab::domain::hmlie M0{y0};

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(5.0);

    lielab::topos::Flow f;
    f.variable_time_step = true;
    f.dt = 0.02;
    lielab::topos::IntegralCurve curve = f(vf, tspan, M0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(std::abs(curve.t(1) - 0.007703769593747) < TOL_COARSE);
    CHECK(std::abs(curve.t(2) - 0.015420629134474) < TOL_COARSE);
    CHECK(std::abs(curve.t(3) - 0.023255332563845) < TOL_COARSE);
    CHECK(std::abs(curve.t(4) - 0.031246577586516) < TOL_COARSE);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.y(0, 0) == y0(0));
    CHECK(curve.y(0, 1) == y0(1));
    CHECK(curve.y(0, 2) == y0(2));

    size_t nrows = curve.y.rows();

    CHECK(std::abs(curve.t(nrows-1) - 5.0) < TOL_FINE);
    CHECK(std::abs(curve.y(nrows-1, 0) - 15.230102737555342) < TOL_COARSE);
    CHECK(std::abs(curve.y(nrows-1, 1) + 0.796697875936802) < TOL_COARSE);
    CHECK(std::abs(curve.y(nrows-1, 2) + 1.472989006310112) < TOL_COARSE);
}

TEST_CASE("Flow_unwrapped_rk45_variable", "[topos]")
{
    /*!
    * Tests Flow against a function that has not been wrapped.
    */

    Eigen::VectorXd y0(3);
    y0(0) = 25.0;
    y0(1) = 0.0;
    y0(2) = -20.0;

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(5.0);

    lielab::topos::Flow f;
    f.variable_time_step = true;
    f.dt = 0.02;
    lielab::topos::IntegralCurve curve = f(eoms_lorenz_unwrapped, tspan, y0);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(std::abs(curve.t(1) - 0.007703769593747) < TOL_COARSE);
    CHECK(std::abs(curve.t(2) - 0.015420629134474) < TOL_COARSE);
    CHECK(std::abs(curve.t(3) - 0.023255332563845) < TOL_COARSE);
    CHECK(std::abs(curve.t(4) - 0.031246577586516) < TOL_COARSE);

    CHECK(std::abs(curve.t(0) - 0.0) < TOL_FINE);
    CHECK(curve.y(0, 0) == y0(0));
    CHECK(curve.y(0, 1) == y0(1));
    CHECK(curve.y(0, 2) == y0(2));

    size_t nrows = curve.y.rows();

    CHECK(std::abs(curve.t(nrows-1) - 5.0) < TOL_FINE);
    CHECK(std::abs(curve.y(nrows-1, 0) - 15.230102737555342) < TOL_COARSE);
    CHECK(std::abs(curve.y(nrows-1, 1) + 0.796697875936802) < TOL_COARSE);
    CHECK(std::abs(curve.y(nrows-1, 2) + 1.472989006310112) < TOL_COARSE);
}

TEST_CASE("A1_vector", "[topos]")
{
    /*!
    * Problem A-1
    *
    * Uses VectorXd representation and events.
    * 
    * Source: Thompson, S. A collection of test problems for ordinary differential
    * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    * Oak Ridge National Lab., TN (USA), 1987.
    */

    lielab::topos::Flow f;

    std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vf = [](const double t, const Eigen::VectorXd & y)
    {
        const double m = 1.0/4.0;
        const double w = 8.0;
        const double k = 2.0;
        Eigen::VectorXd dy = Eigen::VectorXd::Zero(2);
        dy << y(1), (w - k*y(1))/m;
        return dy;
    };

    std::function<double(double, Eigen::VectorXd)> vf_event = [](const double t, const Eigen::VectorXd & y)
    {
        const double H = 10.0;
        return H - y(0);
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(5.0);

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(2);
    y0 << 0.0, 0.0;

    lielab::topos::IntegralCurve out = f(vf, tspan, y0, vf_event);

    const size_t L = out.t.size();

    CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
    CHECK(std::abs(out.t(L-1) - 2.62499999990522) <= 1.0e-9);

    for (size_t ii = 0; ii < L; ii++)
    {
        CHECK(std::abs(out.y(ii, 0) - (4.0*(out.t(ii) + 1.0/8.0*std::exp(-8.0*out.t(ii)) - 1.0/8.0))) <= 1e-6);
        CHECK(std::abs(out.y(ii, 1) - (4.0*(1.0 - std::exp(-8.0*out.t(ii))))) <= 1e-6);
    }
}

TEST_CASE("A1_hom", "[topos]")
{
    /*!
    * Problem A-1
    *
    * Uses hom representation and events.
    * 
    * Source: Thompson, S. A collection of test problems for ordinary differential
    * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    * Oak Ridge National Lab., TN (USA), 1987.
    */

    lielab::topos::Flow f;

    std::function<lielab::domain::halie(double, lielab::domain::hmlie)> vf = [](const double t, const lielab::domain::hmlie & M)
    {
        const double m = 1.0/4.0;
        const double w = 8.0;
        const double k = 2.0;
        const lielab::domain::RN y = std::get<lielab::domain::RN>(M.space[0]);

        const lielab::domain::rn dy{y(1), (w - k*y(1))/m};
        return lielab::domain::halie{dy};
    };

    std::function<double(double, lielab::domain::hmlie)> vf_event = [](const double t, const lielab::domain::hmlie & M)
    {
        const double H = 10.0;
        const lielab::domain::RN y = std::get<lielab::domain::RN>(M.space[0]);
        return H - y(0);
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(5.0);

    const lielab::domain::RN a1y0{0.0, 0.0};
    const lielab::domain::hmlie a1M0{a1y0};

    lielab::topos::IntegralCurve out = f(vf, tspan, a1M0, vf_event);

    const size_t L = out.t.size();

    CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
    CHECK(std::abs(out.t(L-1) - 2.62499999990522) <= 1.0e-9);

    for (size_t ii = 0; ii < L; ii++)
    {
        CHECK(std::abs(out.y(ii, 0) - (4.0*(out.t(ii) + 1.0/8.0*std::exp(-8.0*out.t(ii)) - 1.0/8.0))) <= 1e-6);
        CHECK(std::abs(out.y(ii, 1) - (4.0*(1.0 - std::exp(-8.0*out.t(ii))))) <= 1e-6);
    }
}

TEST_CASE("A2_vector", "[topos]")
{
    /*!
    * Problem A-2
    *
    * Uses VectorXd representation and events.
    * 
    * Source: Thompson, S. A collection of test problems for ordinary differential
    * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    * Oak Ridge National Lab., TN (USA), 1987.
    */

    lielab::topos::Flow f;
    f.dt_max = 10000.0;

    const double g = 32.0;
    const double R = 4000.0 * 5280.0;
    const double H = 237000.0 * 5280.0;

    std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vf = [g, R, H](const double t, const Eigen::VectorXd & y)
    {
        Eigen::VectorXd dy = Eigen::VectorXd::Zero(1);
        dy << -std::sqrt(2.0*g*std::pow(R, 2.0)) * std::sqrt((H - y(0))/(H*y(0)));
        return dy;
    };

    std::function<double(double, Eigen::VectorXd)> vf_event = [g, R, H](const double t, const Eigen::VectorXd & y)
    {
        return y(0) - R;
    };

    std::function<double(double)> h_to_t = [g, R, H](const double h)
    {
        return (std::pow(H, 3.0/2.0) / (8.0*h))*(std::sqrt(h/H - std::pow(h/H, 2.0)) + 1.0/2.0*std::acos(2.0*h/H - 1.0));
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(800000.0);

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    y0 << H - 1.0e-6;
    lielab::topos::IntegralCurve out = f(vf, tspan, y0, vf_event);

    const size_t L = out.t.size();

    CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
    CHECK(std::abs(out.t(L-1) - h_to_t(out.y(L-1, 0))) <= 1.0e-1);
}

TEST_CASE("A2_hom", "[topos]")
{
    /*!
    * Problem A-2
    *
    * Uses hom representation and events.
    * 
    * Source: Thompson, S. A collection of test problems for ordinary differential
    * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    * Oak Ridge National Lab., TN (USA), 1987.
    */

    lielab::topos::Flow f;
    f.dt_max = 10000.0;

    const double g = 32.0;
    const double R = 4000.0 * 5280.0;
    const double H = 237000.0 * 5280.0;

    std::function<lielab::domain::halie(double, lielab::domain::hmlie)> vf = [g, R, H](const double t, const lielab::domain::hmlie & M)
    {
        const lielab::domain::RN y = std::get<lielab::domain::RN>(M.space[0]);
        lielab::domain::rn dy{-std::sqrt(2.0*g*std::pow(R, 2.0)) * std::sqrt((H - y(0))/(H*y(0)))};
        return lielab::domain::halie{dy};
    };

    std::function<double(double, lielab::domain::hmlie)> vf_event = [g, R, H](const double t, const lielab::domain::hmlie & M)
    {
        lielab::domain::RN y = std::get<lielab::domain::RN>(M.space[0]);
        return y(0) - R;
    };

    std::function<double(double)> h_to_t = [g, R, H](const double h)
    {
        return (std::pow(H, 3.0/2.0) / (8.0*h))*(std::sqrt(h/H - std::pow(h/H, 2.0)) + 1.0/2.0*std::acos(2.0*h/H - 1.0));
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(800000.0);

    lielab::domain::RN y0{H - 1.0e-6};
    lielab::topos::IntegralCurve out = f(vf, tspan, lielab::domain::hmlie{y0}, vf_event);

    const size_t L = out.t.size();

    CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
    CHECK(std::abs(out.t(L-1) - h_to_t(out.y(L-1, 0))) <= 1.0e-1);
}

TEST_CASE("A4_vector", "[topos]")
{
    /*!
    * Problem A-4
    *
    * Uses VectorXd representation and events.
    * 
    * Source: Thompson, S. A collection of test problems for ordinary differential
    * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    * Oak Ridge National Lab., TN (USA), 1987.
    */

    lielab::topos::Flow f;
    constexpr double PI = lielab::constants::PI<double>;

    std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vf = [](const double t, const Eigen::VectorXd & y)
    {
        Eigen::VectorXd dy = Eigen::VectorXd::Zero(2);
        dy << y(1), -(16.0 * std::pow(PI, 2) * std::exp(-2.0*t) - 1.0/4.0)*y(0);
        return dy;
    };

    std::function<double(double)> t_to_y1 = [](const double t)
    {
        return std::exp(t/2.0)*std::cos(4.0*PI*std::exp(-t));
    };

    std::function<double(double)> t_to_y2 = [](const double t)
    {
        return std::exp(t/2.0)*(4.0*PI*std::exp(-t)*std::sin(4.0*PI*std::exp(-t)) + 1.0/2.0*std::cos(4.0*PI*std::exp(-t)));
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(std::log(8.0) - std::log(1.0)); // Root at k=1

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(2);
    y0 << 1.0, 0.5;
    
    lielab::topos::IntegralCurve out = f(vf, tspan, y0);

    const size_t L = out.t.size();

    CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
    CHECK(std::abs(out.y(L-1, 0) - 0.0) <= 1.0e-5);

    for (size_t ii = 0; ii < L; ii++)
    {
        CHECK(std::abs(out.y(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-5);
        CHECK(std::abs(out.y(ii, 1) - t_to_y2(out.t(ii))) <= 1.0e-5);
    }
}


TEST_CASE("A4_hom", "[topos]")
{
    /*!
    * Problem A-4
    *
    * Uses hom representation and events.
    * 
    * Source: Thompson, S. A collection of test problems for ordinary differential
    * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    * Oak Ridge National Lab., TN (USA), 1987.
    */

    lielab::topos::Flow f;
    constexpr double PI = lielab::constants::PI<double>;

    std::function<lielab::domain::halie(double, lielab::domain::hmlie)> vf = [](const double t, const lielab::domain::hmlie & M)
    {
        lielab::domain::RN y = std::get<lielab::domain::RN>(M.space[0]);
        lielab::domain::rn dy{y(1), -(16.0 * std::pow(PI, 2) * std::exp(-2.0*t) - 1.0/4.0)*y(0)};
        return lielab::domain::halie{dy};
    };

    std::function<double(double)> t_to_y1 = [](const double t)
    {
        return std::exp(t/2.0)*std::cos(4.0*PI*std::exp(-t));
    };

    std::function<double(double)> t_to_y2 = [](const double t)
    {
        return std::exp(t/2.0)*(4.0*PI*std::exp(-t)*std::sin(4.0*PI*std::exp(-t)) + 1.0/2.0*std::cos(4.0*PI*std::exp(-t)));
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(std::log(8.0) - std::log(1.0)); // Root at k=1

    lielab::domain::RN y0{1.0, 0.5};
    
    lielab::topos::IntegralCurve out = f(vf, tspan, lielab::domain::hmlie{y0});

    const size_t L = out.t.size();

    CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
    CHECK(std::abs(out.y(L-1, 0) - 0.0) <= 1.0e-5);

    for (size_t ii = 0; ii < L; ii++)
    {
        CHECK(std::abs(out.y(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-5);
        CHECK(std::abs(out.y(ii, 1) - t_to_y2(out.t(ii))) <= 1.0e-5);
    }
}

TEST_CASE("A5_vector", "[topos]")
{
    /*!
    * Problem A-5
    *
    * Uses VectorXd representation and events.
    * 
    * Source: Thompson, S. A collection of test problems for ordinary differential
    * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    * Oak Ridge National Lab., TN (USA), 1987.
    */

    lielab::topos::Flow f;
    constexpr double PI = lielab::constants::PI<double>;

    std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vf = [](const double t, const Eigen::VectorXd & y)
    {
        Eigen::VectorXd dy = Eigen::VectorXd::Zero(2);
        dy << y(1), -16.0 * std::cos(PI*t/2.0)*y(1) - (64.0 * std::pow(PI, 2) + 64.0*std::pow(std::cos(PI*t/2.0), 2) - 4.0*PI*std::sin(PI*t/2.0))*y(0);
        return dy;
    };

    std::function<double(double, Eigen::VectorXd)> vf_event = [](const double t, const Eigen::VectorXd & y)
    {
        if ((t > 1.07) && (t < 1.30))
        {
            return -y(0);
        }
        return 1.0;
    };

    std::function<double(double)> t_to_y1 = [](const double t)
    {
        return std::exp(-16/PI*std::sin(PI*t/2))*std::cos(8*PI*t);
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(2.0);

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(2);
    y0 << 1.0, -8.0;
    
    lielab::topos::IntegralCurve out = f(vf, tspan, y0, vf_event);

    const size_t L = out.t.size();

    CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
    CHECK(std::abs(out.t(L-1) - (2.0*10.0-1.0)/16.0) <= 1.0e-5);

    for (size_t ii = 0; ii < L; ii++)
    {
        CHECK(std::abs(out.y(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-5);
    }
}

TEST_CASE("A5_hom", "[topos]")
{
    /*!
    * Problem A-5
    *
    * Uses hom representation and events.
    * 
    * Source: Thompson, S. A collection of test problems for ordinary differential
    * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    * Oak Ridge National Lab., TN (USA), 1987.
    */

    lielab::topos::Flow f;
    constexpr double PI = lielab::constants::PI<double>;

    std::function<lielab::domain::halie(double, lielab::domain::hmlie)> vf = [](const double t, const lielab::domain::hmlie & M)
    {
        lielab::domain::RN y = std::get<lielab::domain::RN>(M.space[0]);
        lielab::domain::rn dy{y(1), -16.0 * std::cos(PI*t/2.0)*y(1) - (64.0 * std::pow(PI, 2) + 64.0*std::pow(std::cos(PI*t/2.0), 2) - 4.0*PI*std::sin(PI*t/2.0))*y(0)};
        return lielab::domain::halie{dy};
    };

    std::function<double(double, lielab::domain::hmlie)> vf_event = [](const double t, const lielab::domain::hmlie & M)
    {
        lielab::domain::RN y = std::get<lielab::domain::RN>(M.space[0]);
        if ((t > 1.07) && (t < 1.30))
        {
            return -y(0);
        }
        return 1.0;
    };

    std::function<double(double)> t_to_y1 = [](const double t)
    {
        return std::exp(-16/PI*std::sin(PI*t/2))*std::cos(8*PI*t);
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(2.0);

    lielab::domain::RN y0{1.0, -8.0};
    
    lielab::topos::IntegralCurve out = f(vf, tspan, lielab::domain::hmlie{y0}, vf_event);

    const size_t L = out.t.size();

    CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-9);
    CHECK(std::abs(out.t(L-1) - (2.0*10.0-1.0)/16.0) <= 1.0e-5);

    for (size_t ii = 0; ii < L; ii++)
    {
        CHECK(std::abs(out.y(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-5);
    }
}

TEST_CASE("B1_vector", "[topos]")
{
    /*!
    * Problem B-1
    *
    * Uses VectorXd representation and events.
    * 
    * Source: Thompson, S. A collection of test problems for ordinary differential
    * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    * Oak Ridge National Lab., TN (USA), 1987.
    */

    lielab::topos::Flow f;

    const double K = 5.0/3.0;
    const double L = 1.0/4.0;

    std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vf = [K, L](const double t, const Eigen::VectorXd & y)
    {
        Eigen::VectorXd dy = Eigen::VectorXd::Zero(2);
        const double dy0 = y(0);
        double dy1 = 0.0;

        const double E = y(0) - y(1);

        if (E < -L/K)
        {
            dy1 = -L;
        }
        else if ((-L/K <= E) && (E <= L/K))
        {
            dy1 = K*E;
        }
        else if (L/K < E)
        {
            dy1 = L;
        }

        dy << dy0, dy1;
        return dy;
    };

    const double t1 = 0.1569;
    const double y2t1 = 1.0199;

    std::function<double(double)> t_to_y1 = [](const double t)
    {
        return std::exp(t);
    };

    std::function<double(double)> t_to_y2 = [K, L, t1, y2t1](const double t)
    {
        if (t < t1)
        {
            return K/(K+1)*std::exp(t) + 1/(K+1)*std::exp(-K*t);
        }
            
        return L*(t - t1) + y2t1;
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(0.5);

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(2);
    y0 << 1.0, 1.0;
    
    lielab::topos::IntegralCurve out = f(vf, tspan, y0);

    const size_t len = out.t.size();

    CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-14);
    CHECK(std::abs(out.t(len-1) - 0.5) <= 1.0e-14);

    for (size_t ii = 0; ii < len; ii++)
    {
        CHECK(std::abs(out.y(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-7);
        CHECK(std::abs(out.y(ii, 1) - t_to_y2(out.t(ii))) <= 1.0e-4); // Answer given in the document is only good to 1e-4 (see t1 and y2(t1))
    }
}

TEST_CASE("B1_hom", "[topos]")
{
    /*!
    * Problem B-1
    *
    * Uses hom representation and events.
    * 
    * Source: Thompson, S. A collection of test problems for ordinary differential
    * equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    * Oak Ridge National Lab., TN (USA), 1987.
    */

    lielab::topos::Flow f;

    std::function<lielab::domain::hmlie(lielab::domain::hmlie, lielab::domain::hmlie)> action = [](const lielab::domain::hmlie & G, const lielab::domain::hmlie & M)
    {
        const lielab::domain::GL _G0 = std::get<lielab::domain::GL>(g.space[0]);
        const lielab::domain::RN _G1 = std::get<lielab::domain::RN>(g.space[1]);
        const lielab::domain::RN _Y0 = std::get<lielab::domain::RN>(M.space[0]);
        const lielab::domain::RN _Y1 = std::get<lielab::domain::RN>(M.space[1]);
        
        const lielab::domain::RN _Y0next = lielab::domain::RN(_G0._data*_Y0._data);
        const lielab::domain::RN _Y1next = _G1*_Y1;
        return lielab::domain::hmlie{_Y0next, _Y1next};
    };

    f.stepper.action = action;

    const double K = 5.0/3.0;
    const double L = 1.0/4.0;

    std::function<lielab::domain::halie(double, lielab::domain::hmlie)> vf = [K, L](const double t, const lielab::domain::hmlie & M)
    {
        lielab::domain::RN _y0 = std::get<lielab::domain::RN>(M.space[0]);
        lielab::domain::RN _y1 = std::get<lielab::domain::RN>(M.space[1]);

        const double y0 = _y0(0);
        const double y1 = _y1(0);

        double dy1 = 0.0;

        const double E = y0 - y1;

        if (E < -L/K)
        {
            dy1 = -L;
        }
        else if ((-L/K <= E) && (E <= L/K))
        {
            dy1 = K*E;
        }
        else if (L/K < E)
        {
            dy1 = L;
        }

        return lielab::domain::halie{lielab::domain::gl::basis(0,1), dy1*lielab::domain::rn::basis(0,2)};
    };

    const double t1 = 0.1569;
    const double y2t1 = 1.0199;

    std::function<double(double)> t_to_y1 = [](const double t)
    {
        return std::exp(t);
    };

    std::function<double(double)> t_to_y2 = [K, L, t1, y2t1](const double t)
    {
        if (t < t1)
        {
            return K/(K+1)*std::exp(t) + 1/(K+1)*std::exp(-K*t);
        }
            
        return L*(t - t1) + y2t1;
    };

    std::vector<double> tspan;
    tspan.push_back(0.0);
    tspan.push_back(0.5);

    lielab::domain::hmlie m0{lielab::domain::RN{1.0}, lielab::domain::RN{1.0}};
    
    lielab::topos::IntegralCurve out = f(vf, tspan, m0);

    const size_t len = out.t.size();

    CHECK(std::abs(out.t(0) - 0.0) <= 1.0e-14);
    CHECK(std::abs(out.t(len-1) - 0.5) <= 1.0e-14);

    for (size_t ii = 0; ii < len; ii++)
    {
        CHECK(std::abs(out.y(ii, 0) - t_to_y1(out.t(ii))) <= 1.0e-14); // Answer is analytic
        CHECK(std::abs(out.y(ii, 1) - t_to_y2(out.t(ii))) <= 1.0e-4); // Answer given in the document is only good to 1e-4 (see t1 and y2(t1))
    }
}
