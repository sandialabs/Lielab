#include <catch2/catch_all.hpp>

#include <Lielab.hpp>

#include <cmath>

TEST_CASE("NewtonMinimize_requires_jacobian", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({2.0, 2.0});

    auto objective = [](const Eigen::VectorXd& x) -> double
    {
        return std::pow(x(0) - 1.0, 2.0) + std::pow(x(1) - 1.0, 2.0);
    };

    EuclideanExtremizationSystem system(objective);

    ExtremizationOptions options;

    NewtonMinimize solver;
    CHECK_THROWS(solver(system, x_guess, options));
}

TEST_CASE("NewtonMinimize_requires_hessian", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({2.0, 2.0});

    auto objective = [](const Eigen::VectorXd& x) -> double
    {
        return std::pow(x[0] - 1.0, 2.0) + std::pow(x[1] - 1.0, 2.0);
    };

    auto jacobian = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd
    {
        return to_VectorXd({2.0*(x[0] - 1.0), 2.0*(x[1] - 1.0)});
    };

    EuclideanExtremizationSystem system(objective);
    system.jacobian = jacobian;

    ExtremizationOptions options;

    NewtonMinimize solver;
    CHECK_THROWS(solver(system, x_guess, options));
}

TEST_CASE("NewtonMinimize_size_hessian", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({2.0, 2.0});

    auto objective = [](const Eigen::VectorXd& x) -> double
    {
        return std::pow(x[0] - 1.0, 2.0) + std::pow(x[1] - 1.0, 2.0);
    };

    auto jacobian = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd
    {
        return to_VectorXd({2.0*(x[0] - 1.0), 2.0*(x[1] - 1.0)});
    };

    auto hessian = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd
    {
        Eigen::MatrixXd out(3,3);
        out(0, 0) = 2.0;
        out(0, 1) = 0.0;
        out(1, 0) = 0.0;
        out(1, 1) = 2.0;
        return out;
    };

    EuclideanExtremizationSystem system(objective);
    system.jacobian = jacobian;
    system.hessian = hessian;

    ExtremizationOptions options;

    NewtonMinimize solver;
    CHECK_THROWS(solver(system, x_guess, options));
}

TEST_CASE("NewtonMinimize_simple", "[optimize]")
{
    using namespace Lielab::optimize;
    using namespace Lielab::utils;

    const Eigen::VectorXd x_guess = to_VectorXd({2.0, 2.0});

    auto objective = [](const Eigen::VectorXd& x) -> double
    {
        return std::pow(x[0] - 1.0, 2.0) + std::pow(x[1] - 1.0, 2.0);
    };

    auto jacobian = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd
    {
        return to_VectorXd({2.0*(x[0] - 1.0), 2.0*(x[1] - 1.0)});
    };

    auto hessian = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd
    {
        Eigen::MatrixXd out(2, 2);
        out(0, 0) = 2.0;
        out(0, 1) = 0.0;
        out(1, 0) = 0.0;
        out(1, 1) = 2.0;
        return out;
    };

    EuclideanExtremizationSystem system(objective);
    system.jacobian = jacobian;
    system.hessian = hessian;

    ExtremizationOptions options;

    NewtonMinimize solver;
    const Eigen::VectorXd xopt = solver(system, x_guess, options);

    REQUIRE(solver.success == true);
    CHECK(std::abs(xopt(0) - 1.0) < options.abstol);
    CHECK(std::abs(xopt(1) - 1.0) < options.abstol);
}

// TEST_CASE("NewtonRootSearch_nans_objective", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.5});

//     auto objective = [](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return to_VectorXd({std::pow(x(0), 2.0)*std::numeric_limits<double>::quiet_NaN() - 1.0});
//     };

//     auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         Eigen::MatrixXd djdx(1, 1);
//         djdx(0, 0) = 2.0*x(0);
//         return djdx;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;

//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == false);
// }

// TEST_CASE("NewtonRootSearch_infs_objective", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.5});

//     auto objective = [](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return to_VectorXd({std::pow(x(0), 2.0)*std::numeric_limits<double>::infinity() - 1.0});
//     };

//     auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         Eigen::MatrixXd djdx(1, 1);
//         djdx(0, 0) = 2.0*x(0);
//         return djdx;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;

//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == false);
// }

// TEST_CASE("NewtonRootSearch_nans_jacobian", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.5});

//     auto objective = [](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
//     };

//     auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         Eigen::MatrixXd djdx(1, 1);
//         djdx(0, 0) = 2.0*x(0)*std::numeric_limits<double>::quiet_NaN();
//         return djdx;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;

//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == false);
// }

// TEST_CASE("NewtonRootSearch_infs_jacobian", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.5});

//     auto objective = [](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
//     };

//     auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         Eigen::MatrixXd djdx(1, 1);
//         djdx(0, 0) = 2.0*x(0)*std::numeric_limits<double>::infinity();
//         return djdx;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;
    
//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == false);
// }

// TEST_CASE("NewtonRootSearch_1D_free", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.5});

//     auto objective = [](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
//     };

//     auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         Eigen::MatrixXd djdx(1, 1);
//         djdx(0, 0) = 2.0*x(0);
//         return djdx;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;

//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == true);
//     CHECK(std::abs(xopt(0) - 1.0) <= options.tol);
// }

// TEST_CASE("NewtonRootSearch_1D_lower", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.5});

//     auto objective = [](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
//     };

//     auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         Eigen::MatrixXd djdx(1, 1);
//         djdx(0, 0) = 2.0*x(0);
//         return djdx;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;
//     system.lower_bound = to_VectorXd({1.1});

//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == false);
//     CHECK(std::abs(xopt(0) - 1.1) <= options.tol);
// }

// TEST_CASE("NewtonRootSearch_1D_upper", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.5});

//     auto objective = [](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
//     };

//     auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         Eigen::MatrixXd djdx(1, 1);
//         djdx(0, 0) = 2.0*x(0);
//         return djdx;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;
//     system.upper_bound = to_VectorXd({0.9});

//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == false);
//     CHECK(std::abs(xopt(0) - 0.9) <= options.tol);
// }

// TEST_CASE("NewtonRootSearch_1D_bounded", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.5});

//     auto objective = [](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return to_VectorXd({std::pow(x(0), 2.0) - 1.0});
//     };

//     auto jacobian = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         Eigen::MatrixXd djdx(1, 1);
//         djdx(0, 0) = 2.0*x(0);
//         return djdx;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;
//     system.lower_bound = to_VectorXd({0.9});
//     system.upper_bound = to_VectorXd({1.1});

//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == true);
//     CHECK(std::abs(xopt(0) - 1.0) <= options.tol);
// }

// TEST_CASE("NewtonRootSearch_2D_free", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.0, 0.0});

//     Eigen::MatrixXd A(2, 2);
//     A(0, 0) = 1.0;
//     A(0, 1) = 2.0;
//     A(1, 0) = -3.0;
//     A(1, 1) = 4.0;

//     Eigen::VectorXd B(2);
//     B(0) = 0.5;
//     B(1) = 0.5;

//     auto objective = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return A*x + B;
//     };

//     auto jacobian = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         return A;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;

//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == true);
//     CHECK(std::abs(xopt(0) + 0.1) <= options.tol);
//     CHECK(std::abs(xopt(1) + 0.2) <= options.tol);
//     CHECK(solver.iteration <= 1);
// }

// TEST_CASE("NewtonRootSearch_2D_lower", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.0, 0.0});

//     Eigen::MatrixXd A(2, 2);
//     A(0, 0) = 1.0;
//     A(0, 1) = 2.0;
//     A(1, 0) = -3.0;
//     A(1, 1) = 4.0;

//     Eigen::VectorXd B(2);
//     B(0) = 0.5;
//     B(1) = 0.5;

//     auto objective = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return A*x + B;
//     };

//     auto jacobian = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         return A;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;
//     system.lower_bound = to_VectorXd({0.0});

//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == false);
//     CHECK(std::abs(xopt(0) + 0.0) <= options.tol);
//     CHECK(std::abs(xopt(1) + 0.2) <= options.tol);
//     CHECK(solver.iteration <= 5);
// }

// TEST_CASE("NewtonRootSearch_3D_underdetermined", "[optimize]")
// {
//     using namespace Lielab::optimize;
//     using namespace Lielab::utils;

//     const Eigen::VectorXd x_guess = to_VectorXd({0.0, 0.0, 0.0});

//     Eigen::MatrixXd A(2, 3);
//     A(0, 0) = 1.0;
//     A(0, 1) = 2.0;
//     A(0, 2) = 4.0;
//     A(1, 0) = -3.0;
//     A(1, 1) = 4.0;
//     A(1, 2) = -5.0;

//     Eigen::VectorXd B(2);
//     B(0) = 0.5;
//     B(1) = 0.5;

//     auto objective = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd
//     {
//         return A*x + B;
//     };

//     auto jacobian = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd
//     {
//         return A;
//     };

//     EuclideanRootSystem system(objective);
//     system.jacobian = jacobian;

//     RootSolvingOptions options;

//     NewtonRootSearch solver;
//     const Eigen::VectorXd xopt = solver(system, x_guess, options);

//     CHECK(solver.success == true);
//     const Eigen::VectorXd objfinal = objective(xopt);
//     CHECK(std::abs(objfinal(0)) <= options.tol);
//     CHECK(std::abs(objfinal(1)) <= options.tol);
//     CHECK(solver.iteration <= 1);
// }
