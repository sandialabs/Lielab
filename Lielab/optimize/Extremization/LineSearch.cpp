#include "LineSearch.hpp"

#include "ExtremizationCommon.hpp"

#include <Eigen/Core>

namespace Lielab::optimize
{

LineSearch::LineSearch()
{

}

Eigen::VectorXd LineSearch::operator()(const EuclideanExtremizationSystem system,
                                       const Eigen::VectorXd& guess,
                                       const Eigen::VectorXd& direction,
                                       const double magnitude,
                                       const ExtremizationOptions options)
{
    this->iteration = 0;
    // this->status = 
    this->success = false;
    this->message = "";

    const int max_iter = (options.max_iterations < 0) ? 10 : options.max_iterations;

    if (guess.size() != direction.size())
    {
        // this->status = 
        this->success = false;
        this->message = "LineSearch: Size of guess must be equal to size of direction.";
        return guess;
    }

    const double dx = 1.0e-10;

    // TODO: Make these not numerical?
    const double f0 = system.objective(guess);
    const double f1 = system.objective(guess + dx*direction);
    double dfdx = (f1 - f0)/dx;

    Eigen::VectorXd _direction = direction;
    if (dfdx > 0)
    {
        // Reverse the direction if dfdx is increasing
        _direction.noalias() = -direction;
        dfdx = -dfdx;
    }

    // Initialize alpha
    this->alpha = options.initial_step_size/_direction.norm();
    if (!std::isnan(options.initial_alpha))
    {
        this->alpha = options.initial_alpha;
    }

    Eigen::VectorXd xnew;
    double fnew;
    Eigen::VectorXd xbest = guess;
    double fbest = f0;

    while (this->iteration < max_iter)
    {
        this->iteration++;

        // Create new guess
        xnew.noalias() = guess + this->alpha*_direction;

        // Account for bounds
        for (int ii = 0; ii < static_cast<int>(system.lower_bound.size()); ii++)
        {
            if (!std::isnan(system.lower_bound(ii)))
            {
                xnew(ii) = std::max(xnew(ii), system.lower_bound(ii));
            }
        }

        for (int ii = 0; ii < static_cast<int>(system.upper_bound.size()); ii++)
        {
            if (!std::isnan(system.upper_bound(ii)))
            {
                xnew(ii) = std::min(xnew(ii), system.upper_bound(ii));
            }
        }

        fnew = system.objective(xnew);

        // if (fnew <= f0 - options.sufficient_decrease*this->alpha*std::abs(magnitude)) // TODO: Use this factor
        if (fnew <= f0 + options.sufficient_decrease*this->alpha*dfdx)
        {
            this->success = true;
            this->message = "LineSearch: Converged to tolerance.";
            return xnew;
        }
        else if (fnew < fbest)
        {
            fbest = fnew;
            xbest.noalias() = xnew;
        }

        this->alpha *= options.contraction_factor;
    }

    this->success = false;
    this->message = "LineSearch: Max iterations exceeded.";
    return xbest;
}

}
