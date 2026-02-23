#include <catch2/catch_all.hpp>

#include <Lielab/optimize.hpp>
#include <Lielab/utils.hpp>

#include <cmath>

TEST_CASE("HybridRootSearch_requires_jacobian", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({0.5});

    auto objective = [](const Eigen::VectorXd& x) -> Eigen::VectorXd
    {
        return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
    };

    EuclideanRootSystem system(objective);

    RootSolvingOptions options;

    HybridRootSearch solver;
    const Eigen::VectorXd xopt = solver(system, x_guess, options);

    CHECK(solver.success == false);
}

TEST_CASE("HybridRootSearch_no_bounds", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({0.5});

    auto objective = [](const Eigen::VectorXd& x)
    {
        return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
    };

    auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
    {
        Eigen::MatrixXd djdx(1, 1);
        djdx(0, 0) = 2.0*x(0);
        return djdx;
    };

    EuclideanRootSystem system(objective);
    system.jacobian = jacobian;
    system.lower_bound = to_VectorXd({-2.0});
    system.upper_bound = to_VectorXd({2.0});

    RootSolvingOptions options;

    HybridRootSearch solver;
    const Eigen::VectorXd xopt = solver(system, x_guess, options);

    CHECK(solver.success == false);
}

TEST_CASE("HybridRootSearch_nans_objective", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({0.5});

    auto objective = [](const Eigen::VectorXd& x)
    {
        return to_VectorXd({std::pow(x(0), 2.0)*std::numeric_limits<double>::quiet_NaN() - 1.0});
    };

    auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
    {
        Eigen::MatrixXd djdx(1, 1);
        djdx(0, 0) = 2.0*x(0);
        return djdx;
    };

    EuclideanRootSystem system(objective);
    system.jacobian = jacobian;

    RootSolvingOptions options;

    HybridRootSearch solver;
    const Eigen::VectorXd xopt = solver(system, x_guess, options);

    CHECK(solver.success == false);
}

TEST_CASE("HybridRootSearch_infs_objective", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({0.5});

    auto objective = [](const Eigen::VectorXd& x)
    {
        return to_VectorXd({std::pow(x(0), 2.0)*std::numeric_limits<double>::infinity() - 1.0});
    };

    auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
    {
        Eigen::MatrixXd djdx(1, 1);
        djdx(0, 0) = 2.0*x(0);
        return djdx;
    };

    EuclideanRootSystem system(objective);
    system.jacobian = jacobian;

    RootSolvingOptions options;

    HybridRootSearch solver;
    const Eigen::VectorXd xopt = solver(system, x_guess, options);

    CHECK(solver.success == false);
}

TEST_CASE("HybridRootSearch_nans_jacobian", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({0.5});

    auto objective = [](const Eigen::VectorXd& x)
    {
        return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
    };

    auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
    {
        Eigen::MatrixXd djdx(1, 1);
        djdx(0, 0) = 2.0*x(0)*std::numeric_limits<double>::quiet_NaN();
        return djdx;
    };

    EuclideanRootSystem system(objective);
    system.jacobian = jacobian;

    RootSolvingOptions options;

    HybridRootSearch solver;
    const Eigen::VectorXd xopt = solver(system, x_guess, options);

    CHECK(solver.success == false);
}

TEST_CASE("HybridRootSearch_infs_jacobian", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({0.5});

    auto objective = [](const Eigen::VectorXd& x)
    {
        return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
    };

    auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
    {
        Eigen::MatrixXd djdx(1, 1);
        djdx(0, 0) = 2.0*x(0)*std::numeric_limits<double>::infinity();
        return djdx;
    };

    EuclideanRootSystem system(objective);
    system.jacobian = jacobian;

    RootSolvingOptions options;

    HybridRootSearch solver;
    const Eigen::VectorXd xopt = solver(system, x_guess, options);

    CHECK(solver.success == false);
}

TEST_CASE("HybridRootSearch_1D_free", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({0.5});

    auto objective = [](const Eigen::VectorXd& x)
    {
        return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
    };

    auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
    {
        Eigen::MatrixXd djdx(1, 1);
        djdx(0, 0) = 2.0*x(0);
        return djdx;
    };

    EuclideanRootSystem system(objective);
    system.jacobian = jacobian;

    RootSolvingOptions options;

    HybridRootSearch solver;
    const Eigen::VectorXd xopt = solver(system, x_guess, options);

    CHECK(solver.success == true);
    CHECK(std::abs(xopt(0) - 1.0) <= options.tol);
}

TEST_CASE("HybridRootSearch_2D", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({0.0, 0.0});

    Eigen::MatrixXd A(2, 2);
    A(0, 0) = 1.0;
    A(0, 1) = 2.0;
    A(1, 0) = -3.0;
    A(1, 1) = 4.0;

    Eigen::VectorXd B(2);
    B(0) = 0.5;
    B(1) = 0.5;

    auto objective = [&](const Eigen::VectorXd& x)
    {
        return A*x + B;
    };

    auto jacobian = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd
    {
        return A;
    };

    EuclideanRootSystem system(objective);
    system.jacobian = jacobian;

    RootSolvingOptions options;

    HybridRootSearch solver;
    const Eigen::VectorXd xopt = solver(system, x_guess, options);

    CHECK(solver.success == true);
    CHECK(std::abs(xopt(0) + 0.1) <= options.tol);
    CHECK(std::abs(xopt(1) + 0.2) <= options.tol);
}

TEST_CASE("HybridRootSearch_3D_underdetermined", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({0.0, 0.0, 0.0});

    Eigen::MatrixXd A(2, 3);
    A(0, 0) = 1.0;
    A(0, 1) = 2.0;
    A(0, 2) = 4.0;
    A(1, 0) = -3.0;
    A(1, 1) = 4.0;
    A(1, 2) = -5.0;

    Eigen::VectorXd B(2);
    B(0) = 0.5;
    B(1) = 0.5;

    auto objective = [&](const Eigen::VectorXd& x)
    {
        return A*x + B;
    };

    auto jacobian = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd
    {
        return A;
    };

    EuclideanRootSystem system(objective);
    system.jacobian = jacobian;

    RootSolvingOptions options;

    HybridRootSearch solver;
    const Eigen::VectorXd xopt = solver(system, x_guess, options);

    CHECK(solver.success == false);
    // const Eigen::VectorXd objfinal = objective(xopt.x);
    // CHECK(std::abs(objfinal(0)) <= options.tol);
    // CHECK(std::abs(objfinal(1)) <= options.tol);
}
