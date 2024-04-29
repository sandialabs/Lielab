
double optimfun(Eigen::VectorXd x)
{
    double out = (x(0)*x(0)*x(0) - 2*x(0) - 5);
    return out;
}

TEST_CASE("opt_golden", "[optim]")
{
    Lielab::optim::opt_golden search;
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

    Eigen::VectorXd opt = search(*optimfun);

    CHECK(std::abs(opt(0) - 0.816496623231839) < TOL_FINE);
    CHECK(search.iterations == 31);
    CHECK(search.success == true);
    CHECK(search.algo_status == Lielab::ALGO_STATUS::FINISHED);

    search.init();
    search.max_iterations = 5;

    opt = search(*optimfun);

    CHECK(search.success == false);
    CHECK(search.iterations == 5);
    CHECK(search.algo_status == Lielab::ALGO_STATUS::MAXITER);

}

double optfun2(const Lielab::domain::CompositeAlgebra & m)
{
    Lielab::domain::rn y = std::get<Lielab::domain::rn>(m.space[0]);
    return std::pow(y(0), 2) - 1;
}

TEST_CASE("hnewton", "[optim]")
{
    Eigen::VectorXd _y01(1);
    _y01 << 0.5;
    Lielab::domain::rn _y02(_y01);
    Lielab::domain::CompositeAlgebra y01{_y02};

    Lielab::optim::hnewton search;

    Lielab::domain::CompositeAlgebra yopt = search(optfun2, y01);
    Lielab::domain::rn yopt2 = std::get<Lielab::domain::rn>(yopt.space[0]);

    CHECK(std::abs(yopt2(0) - 1.0) <= TOL_FINE);
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
