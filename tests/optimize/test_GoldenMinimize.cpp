#include <catch2/catch_all.hpp>
#include <Eigen/Core>

#include <Lielab.hpp>
#include "../test_utils.hpp"

TEST_CASE("GoldenMinimize_fail_bounds", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    auto optimfun = [&](const Eigen::VectorXd& x)
    {
        return (x(0)*x(0)*x(0) - 2*x(0) - 5);
    };

    EuclideanExtremizationSystem system(optimfun);

    ExtremizationOptions options;

    GoldenMinimize solver;
    const Eigen::VectorXd xopt = solver(system, options);

    CHECK(solver.success == false);
}

TEST_CASE("GoldenMinimize_1D", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    auto optimfun = [&](const Eigen::VectorXd& x)
    {
        return (x(0)*x(0)*x(0) - 2*x(0) - 5);
    };

    EuclideanExtremizationSystem system(optimfun);
    system.lower_bound = to_VectorXd({0.0});
    system.upper_bound = to_VectorXd({2.0});

    ExtremizationOptions options;

    GoldenMinimize solver;
    const Eigen::VectorXd xopt = solver(system, options);

    CHECK(solver.success == true);
    CHECK(std::abs(xopt(0) - std::sqrt(2.0/3.0)) < options.reltol);
}

TEST_CASE("GoldenMinimize_1D_right_nans", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    auto optimfun = [&](const Eigen::VectorXd& x)
    {
        if (x(0) > 0.8) return std::numeric_limits<double>::quiet_NaN();
        return 1.0 - x(0);
    };

    EuclideanExtremizationSystem system(optimfun);
    system.lower_bound = to_VectorXd({0.0});
    system.upper_bound = to_VectorXd({2.0});

    ExtremizationOptions options;

    GoldenMinimize solver;
    const Eigen::VectorXd xopt = solver(system, options);

    CHECK(solver.success == true);
    CHECK(std::abs(xopt(0) - 0.8) < options.reltol);
}

TEST_CASE("GoldenMinimize_1D_left_nans", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    auto optimfun = [&](const Eigen::VectorXd& x)
    {
        if (x(0) < 1.2) return std::numeric_limits<double>::quiet_NaN();
        return x(0);
    };

    EuclideanExtremizationSystem system(optimfun);
    system.lower_bound = to_VectorXd({0.0});
    system.upper_bound = to_VectorXd({2.0});

    ExtremizationOptions options;

    GoldenMinimize solver;
    const Eigen::VectorXd xopt = solver(system, options);

    CHECK(solver.success == true);
    CHECK(std::abs(xopt(0) - 1.2) < options.reltol);
}
