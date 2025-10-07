#ifndef LIELAB_UTILS_NEWTON_HPP
#define LIELAB_UTILS_NEWTON_HPP

#include "Error.hpp"

#include <Eigen/Core>

namespace Lielab::utils
{

class newton
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

    // Newton vals
    double fdt = 1e-6;
    double lower;
    double upper;

    double f;
    double fnext;
    double m;
    double mnext;

    size_t _dim;

    newton();
    int init(const double m0, double fval);
    int step0();
    int step1();
    double operator()(const std::function<double(double)> & fun, const double m0);

};

}

#endif
