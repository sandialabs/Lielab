#include "GoldenMinimize.hpp"

#include "ExtremizationCommon.hpp"

#include <Eigen/Core>

namespace Lielab::optimize
{

GoldenMinimize::GoldenMinimize()
{

}

Eigen::VectorXd GoldenMinimize::operator()(EuclideanExtremizationSystem system, const ExtremizationOptions options)
{
    this->iteration = 0;
    this->success = false;
    this->message = "";

    if (system.lower_bound.size() != system.upper_bound.size())
    {
        this->message = "GoldenMinimize: lower_bound and upper_bound must be equal size.";
        this->success = false;
        return system.lower_bound;
    }

    if (system.lower_bound.size() != 1)
    {
        this->message = "GoldenMinimize: Only 1D systems are implemented.";
        this->success = false;
        return system.lower_bound;
    }

    if (system.lower_bound.array().isNaN().any() || system.upper_bound.array().isNaN().any())
    {
        this->message = "GoldenMinimize: NaNs found in lower_bound or upper_bound.";
        this->success = false;
        return system.lower_bound;
    }

    if (system.lower_bound.array().isInf().any() || system.upper_bound.array().isInf().any())
    {
        this->message = "GoldenMinimize: Infs found in lower_bound or upper_bound.";
        this->success = false;
        return system.lower_bound;
    }

    const double tau = 1.0/options.golden_ratio;

    Eigen::VectorXd x0 = system.lower_bound;
    Eigen::VectorXd x3 = system.upper_bound;

    Eigen::VectorXd x1;
    Eigen::VectorXd x2;

    // TODO: The start could be made more robust to nans if
    //       we do a manual check and start here.
    // Eigen::VectorXd x1 = (x0 + x3)/2.0;
    // Eigen::VectorXd x2 = x1 + (1.0-tau)*(x3 - x1);
    x1 = x0 + (1.0 - tau)*(x3 - x0);
    x2 = x0 + tau*(x3 - x0);
    
    // const double dist0 = (xright - xleft).norm();

    double f1 = system.objective(x1);
    double f2 = system.objective(x2);

    if ((std::isnan(f1) || std::isinf(f1)) && (std::isnan(f2) || std::isinf(f2)))
    {
        this->message = "GoldenMinimize: Too many Infs or NaNs in objective function.";
        return system.lower_bound;
    }

    bool running = true;
    if (options.max_iterations == 0) running = false;

    while (running)
    {
        this->iteration++;
        
        if ((std::isnan(f2) || std::isinf(f2)) || f1 < f2)
        {
            x3 = x2;
            x2 = x1;
            x1 = tau*x2 + (1.0-tau)*x0;

            f2 = f1;
            f1 = system.objective(x1);
            // xright = x2;
            // x2 = x1;
            // f2 = f1;
            // x1 = xleft + (1.0 - this->tau)*(xright - xleft);
            // f1 = objective(x1);
        }
        else // Check isinf/nan of f1?
        {
            x0 = x1;
            x1 = x2;
            x2 = tau*x1 + (1.0-tau)*x3;

            f1 = f2;
            f2 = system.objective(x2);
            // xleft = x1;
            // x1 = x2;
            // f1 = f2;
            // x2 = xleft + this->tau*(xright - xleft);
            // f2 = objective(x2);
        }

        if ((x2 - x1).norm() < options.abstol)
        {
            this->success = true;
            this->message = "GoldenMinimize: Converged to absolute tolerance.";
            running = false;
        }

        if ((x3 - x0).norm() < options.reltol*(x1.norm() + x3.norm()))
        {
            this->success = true;
            this->message = "GoldenMinimize: Converged to relative tolerance.";
            running = false;
        }

        if (options.max_iterations > 0 && this->iteration >= options.max_iterations)
        {
            this->success = true;
            this->message = "GoldenMinimize: Converged to relative tolerance.";
            running = false;
        }
    }

    const double fmid = system.objective((x1 + x2)/2.0);

    this->success = true;
    if (f1 < f2 && f1 < fmid)
    {
        return x1;
    }
    else if (f2 < f1 && f2 < fmid)
    {
        return x2;
    }

    return (x1 + x2)/2.0;
}

}
