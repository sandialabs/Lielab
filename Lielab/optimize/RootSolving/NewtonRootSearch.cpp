#include "NewtonRootSearch.hpp"
#include "RootSolvingCommon.hpp"

#include <Lielab/utils.hpp>

#include <Eigen/Dense>
#include <Eigen/QR>

#include <functional>
#include <limits>

#include <iostream>

namespace Lielab::optimize
{

NewtonRootSearch::NewtonRootSearch()
{

}

Eigen::VectorXd NewtonRootSearch::operator()(const EuclideanRootSystem system, const Eigen::VectorXd& x_guess, const RootSolvingOptions options)
{
    RootSolvingOptions _options = options;

    Eigen::VectorXd xcurrent = x_guess;
    Eigen::VectorXd fcurrent = system.objective(xcurrent);

    // TODO: Check for infs and nans in fcurrent here?

    Eigen::MatrixXd dfdx = system.jacobian(x_guess);
    if (dfdx.rows() == 1 && dfdx.cols() == 2 && std::isnan(dfdx(0, 0)) && std::isnan(dfdx(0, 1)))
    {
        this->message = "NewtonRootSearch: Jacobian function must be defined.";
        this->success = false;
        return xcurrent;
    }

    const ptrdiff_t n_vars = xcurrent.size();
    const ptrdiff_t n_funcs = fcurrent.size();
    const ptrdiff_t n_lower = system.lower_bound.size();
    const ptrdiff_t n_upper = system.upper_bound.size();

    if (n_vars < n_funcs)
    {
        // Overdetermined systems not implemented.
        this->message = "NewtonRootSearch: Number of free variables must be equal to or greater than functions.";
        this->success = false;
        return xcurrent;
    }

    auto _obj = [&](const Eigen::VectorXd& x)
    {
        this->num_objective_evals++;
        return system.objective(x);
    };
    
    Eigen::VectorXd xnext(n_vars);
    Eigen::VectorXd xdelta(n_vars);
    Eigen::VectorXd diff(n_vars);
    
    this->iteration = 0;
    this->success = false;
    this->message = "";

    // Fix initial guess if it's outside bounds.
    for (ptrdiff_t ii = 0; ii < n_lower; ii++)
    {
        if (!std::isnan(system.lower_bound(ii)) && xcurrent(ii) < system.lower_bound(ii))
        {
            xcurrent(ii) = system.lower_bound(ii);
        }
    }

    for (ptrdiff_t ii = 0; ii < n_upper; ii++)
    {
        if (!std::isnan(system.upper_bound(ii)) && xcurrent(ii) > system.upper_bound(ii))
        {
            xcurrent(ii) = system.upper_bound(ii);
        }
    }

    while (this->success == false)
    {
        // Pre-iteration termination conditions
        if ((_options.max_iterations != -999) && iteration >= _options.max_iterations)
        {
            // End the process if we exceed max iterations. Ignored when max_iterations is -999.
            this->message = "NewtonRootSearch: Max iterations reached.";
            this->success = false;
            return xcurrent;
        }

        fcurrent = _obj(xcurrent);

        if (fcurrent.array().isNaN().any())
        {
            this->message = "NewtonRootSearch: NaNs in objective.";
            this->success = false;
            return x_guess;
        }

        if (fcurrent.array().isInf().any())
        {
            this->message = "NewtonRootSearch: Infs in objective.";
            this->success = false;
            return x_guess;
        }

        if ((fcurrent.cwiseAbs().array() < _options.tol).all()) // TODO: Or there are more zeros than variables
        {
            // End if the current x value satisfies f(x) = 0.
            this->message = "NewtonRootSearch: Converged to tolerance.";
            this->success = true;
            return xcurrent;
        }

        // Do the Newton iteration.
        this->iteration++;
        dfdx = system.jacobian(xcurrent);

        if (dfdx.array().isNaN().any())
        {
            this->message = "NewtonRootSearch: NaNs in jacobian.";
            this->success = false;
            return x_guess;
        }

        if (dfdx.array().isInf().any())
        {
            this->message = "NewtonRootSearch: Infs in jacobian.";
            this->success = false;
            return x_guess;
        }
        
        if (n_vars == n_funcs)
        {
            // System is well determined.
            xdelta = dfdx.inverse()*(-fcurrent);
            // TODO: Check if the inverse fails, if it does try pinv?
        }
        else if (n_vars > n_funcs)
        {
            // System is under determined.
            xdelta = dfdx.completeOrthogonalDecomposition().pseudoInverse()*(-fcurrent);
        }
        else
        {
            // System is over determined
            // xdelta = xcurrent;
        }

        // TODO: Check if constraints are active and do a smarter step.
        
        xnext = xcurrent + xdelta;

        for (ptrdiff_t ii = 0; ii < n_lower; ii++)
        {
            if (!std::isnan(system.lower_bound(ii)) && xnext(ii) < system.lower_bound(ii))
            {
                xnext(ii) = system.lower_bound(ii);
            }
        }

        for (ptrdiff_t ii = 0; ii < n_upper; ii++)
        {
            if (!std::isnan(system.upper_bound(ii)) && xnext(ii) > system.upper_bound(ii))
            {
                xnext(ii) = system.upper_bound(ii);
            }
        }

        diff = xnext - xcurrent;
        
        if ((diff.cwiseAbs().array() < _options.tol*0.1).all())
        {
            // End if progress is not being made.
            this->message = "NewtonRootSearch: Progress is not being made (change in x too small).";
            this->success = false;
            return xnext;
        }

        xcurrent = xnext;
    }

    // This shouldnt be called.
    return Lielab::utils::to_VectorXd({std::numeric_limits<double>::quiet_NaN()});
}

}
