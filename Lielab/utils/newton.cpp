#include "newton.hpp"

#include "Error.hpp"

#include <Eigen/Core>

namespace Lielab::utils
{

newton::newton()
{

}

int newton::init(const double m0, double fval)
{
    this->iterations = 0;
    this->success = false;
    
    this->num_objective_evals = 0;
    this->num_jacobian_evals = 0;
    this->num_hessian_evals = 0;

    this->fnext = fval;
    this->mnext = m0;
    return 2;
}

int newton::step0()
{
    this->iterations++;

    if (max_iterations > 0 && iterations >= max_iterations)
    {
        return 1;
    }

    this->f = this->fnext;
    this->m = this->mnext;

    double dm = this->fdt;

    this->mnext = this->m + dm;

    return 2;
}

int newton::step1()
{
    const double fprime = (this->fnext - this->f)/(this->fdt);

    double dm = -this->f / fprime;

    this->mnext = this->m + dm;

    dm = this->mnext - this->m;
    double err = dm;

    if (std::abs(err) < tolerance)
    {
        success = true;
        return 0;
    }

    return 2;
}

double newton::operator()(const std::function<double(double)> & fun, const double m0)
{
    int status = init(m0, fun(m0));

    while (status != 0 && status != 1)
    {
        status = step0();
        this->fnext = fun(this->mnext);
        status = step1();
        this->fnext = fun(this->mnext);
    }

    this->m = this->mnext;

    return this->m;
}

}
