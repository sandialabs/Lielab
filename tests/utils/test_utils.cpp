#include <catch2/catch_all.hpp>
#include <Eigen/Core>

#include <Lielab.hpp>
#include "../test_utils.hpp"

double optimfun(Eigen::VectorXd x)
{
    double out = (x(0)*x(0)*x(0) - 2*x(0) - 5);
    return out;
}

TEST_CASE("opt_golden", "[utils]")
{
    Lielab::utils::opt_golden search;
    Eigen::VectorXd lower(1);
    lower(0) = 0.0;
    Eigen::VectorXd upper(1);
    upper(0) = 2.0;

    search.lower = lower;
    search.upper = upper;

    search.init();

    CHECK(std::abs(search.lower(0) - 0.0) < TOL_FINE);
    CHECK(std::abs(search.upper(0) - 2.0) < TOL_FINE);
    CHECK(std::abs(search._A(0) - 0.0) < TOL_FINE);
    CHECK(std::abs(search._B(0) - 2.0) < TOL_FINE);
    CHECK(std::abs(search._X1(0) - 0.763932022500210) < TOL_FINE);
    CHECK(std::abs(search._X2(0) - 1.236067977499790) < TOL_FINE);

    search._f1 = optimfun(search._X1);
    search.num_objective_evals++;
    search._f2 = optimfun(search._X2);
    search.num_objective_evals++;
    search.step();

    CHECK(std::abs(search.lower(0) - 0.0) < TOL_FINE);
    CHECK(std::abs(search.upper(0) - 2.0) < TOL_FINE);
    CHECK(std::abs(search._A(0) - 0.0) < TOL_FINE);
    CHECK(std::abs(search._B(0) - 1.236067977499790) < TOL_FINE);
    CHECK(std::abs(search._X1(0) - 0.472135954999579) < TOL_FINE);
    CHECK(std::abs(search._X2(0) - 0.763932022500210) < TOL_FINE);

    search._f1 = optimfun(search._X1);
    search.num_objective_evals++;
    search._f2 = optimfun(search._X2);
    search.num_objective_evals++;
    search.step();

    CHECK(std::abs(search.lower(0) - 0.0) < TOL_FINE);
    CHECK(std::abs(search.upper(0) - 2.0) < TOL_FINE);
    CHECK(std::abs(search._A(0) - 0.472135954999579) < TOL_FINE);
    CHECK(std::abs(search._B(0) - 1.236067977499790) < TOL_FINE);
    CHECK(std::abs(search._X1(0) - 0.763932022500210) < TOL_FINE);
    CHECK(std::abs(search._X2(0) - 0.944271909999159) < TOL_FINE);

    search.init();

    CHECK(search.iterations == 0);
    CHECK(search.num_objective_evals == 0);
    CHECK(search.num_jacobian_evals == 0);
    CHECK(search.num_hessian_evals == 0);

    CHECK(std::abs(search.lower(0) - 0.0) < TOL_FINE);
    CHECK(std::abs(search.upper(0) - 2.0) < TOL_FINE);
    CHECK(std::abs(search._A(0) - 0.0) < TOL_FINE);
    CHECK(std::abs(search._B(0) - 2.0) < TOL_FINE);
    CHECK(std::abs(search._X1(0) - 0.763932022500210) < TOL_FINE);
    CHECK(std::abs(search._X2(0) - 1.236067977499790) < TOL_FINE);

    Lielab::utils::golden_solution opt = search(*optimfun);

    CHECK(std::abs(opt.x(0) - 0.816496623231839) < TOL_FINE);
    CHECK(opt.iterations == 31);
    CHECK(opt.status == 0);

    search.init();
    search.max_iterations = 5;

    opt = search(*optimfun);

    CHECK(opt.iterations == 5);
    CHECK(opt.status == 1);

}

double optfun2(const double x)
{
    return std::pow(x, 2.0) - 1;
}

TEST_CASE("newton", "[utils]")
{
    double y01 = 0.5;

    Lielab::utils::newton search;

    double yopt = search(optfun2, y01);

    CHECK(std::abs(yopt - 1.0) <= TOL_FINE);
}

double myfun(double x)
{
    return x - 4.5;
}

// TEST_CASE("search_linearx", "[optim]")
// {
//     Lielab::optim::search_linearx search;
//     double x = 5.0;

//     search.lower = 4.0;
//     search.upper = 6.0;

//     search.init(x);

//     CHECK(std::abs(search.lower - 4.0) < TOL_FINE);
//     CHECK(std::abs(search.upper - 6.0) < TOL_FINE);

//     x = search.step(x, myfun(x));

//     CHECK(std::abs(x - 5.000006) < TOL_FINE);
//     CHECK(search.k == 3);

//     x = search.step(x, myfun(x));

//     CHECK(std::abs(x - 4.5) < TOL_FINE);
//     CHECK(search.k == 4);

//     x = search.step(x, myfun(x));

//     CHECK(std::abs(x - 4.5) < TOL_FINE);
//     CHECK(search.k == 8);

//     x = 5.0;

//     x = search(*myfun, x);

//     CHECK(std::abs(x - 4.5) < TOL_FINE);
//     CHECK(search.k == 8);

// }
