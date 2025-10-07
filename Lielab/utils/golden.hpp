#ifndef LIELAB_UTILS_GOLDEN_HPP
#define LIELAB_UTILS_GOLDEN_HPP

#include <Eigen/Core>

#include <cmath>

namespace Lielab::utils
{

struct golden_solution
{
    int iterations;
    int status;
    double fval;
    Eigen::VectorXd x;
};

class opt_golden
{
    public:

    // general
    int iterations = 0;
    int max_iterations = -1;
    bool success = false;

    // nlp specific

    int num_objective_evals = 0;
    int num_jacobian_evals = 0;
    int num_hessian_evals = 0;
    double tolerance = 1e-6;

    double val_objective;
    Eigen::VectorXd val_jacobian;
    Eigen::MatrixXd val_hessian;
    
    Eigen::VectorXd x0;
    Eigen::VectorXd lower;
    Eigen::VectorXd upper;

    // Golden specific

    double tau = (std::sqrt(5.0)-1.0)/2.0;

    double _f1;
    double _f2;
    Eigen::VectorXd _X;
    Eigen::VectorXd _X1;
    Eigen::VectorXd _X2;
    Eigen::VectorXd _A;
    Eigen::VectorXd _B;

    opt_golden();
    int init();
    int step();
    golden_solution operator()(double (*obj_fun)(Eigen::VectorXd));

};

}

#endif
