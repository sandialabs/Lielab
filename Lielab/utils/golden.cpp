#include "golden.hpp"

#include <Eigen/Core>

#include <cmath>

namespace Lielab::utils
{

opt_golden::opt_golden()
{

}

int opt_golden::init()
{
    this->iterations = 0;
    this->success = false;
    
    this->num_objective_evals = 0;
    this->num_jacobian_evals = 0;
    this->num_hessian_evals = 0;

    _X1 = lower + (1-tau)*(upper - lower);
    _X2 = lower + tau*(upper - lower);
    _A = lower;
    _B = upper;

    return 2;
}

int opt_golden::step()
{
    /*!
    * Return flags:
    *     - 0: Optimization finished successfully.
    *     - 1: Optimization finished, max iterations exceeded.
    *     - 2: Optimization running, next iteration.
    *
    */

    this->iterations++;

    if (max_iterations > 0 && iterations >= max_iterations)
    {
        return 1;
    }

    if (this->_f1 < this->_f2)
    {
        this->val_objective = this->_f1;
    }
    else
    {
        this->val_objective = this->_f2;
    }
    
    if (this->_f1 < this->_f2)
    {
        this->_B = this->_X2;
        this->_X2 = this->_X1;
        this->_X1 = this->_A + (1.0 - this->tau)*(this->_B - this->_A);
    }
    else
    {
        this->_A = this->_X1;
        this->_X1 = this->_X2;
        this->_X2 = this->_A + this->tau*(this->_B - this->_A);
    }

    if (std::abs((this->_B - this->_A).norm()) < this->tolerance)
    {
        this->success = true;
        return 0;
    }

    return 2;
}

golden_solution opt_golden::operator()(double (*obj_fun)(Eigen::VectorXd))
{
    int status = this->init();

    while (status != 0 && status != 1)
    {
        this->_f1 = obj_fun(this->_X1);
        this->num_objective_evals++;

        this->_f2 = obj_fun(this->_X2);
        this->num_objective_evals++;
        status = this->step();
    }
    
    this->_f1 = obj_fun(_X1);
    this->num_objective_evals++;

    this->_f2 = obj_fun(_X2);
    this->num_objective_evals++;

    if (this->_f1 < this->_f2)
    {
        this->val_objective = this->_f1;
        this->_X = this->_X1;
    }
    else
    {
        this->val_objective = this->_f2;
        this->_X = this->_X2;
    }

    golden_solution out;
    out.iterations = this->iterations;
    out.fval = this->val_objective;
    out.status = status;
    out.x = this->_X;

    return out;
}

}
