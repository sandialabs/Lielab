#include "NewtonMinimize.hpp"

#include "ExtremizationCommon.hpp"
#include "LineSearch.hpp"

#include <Lielab/testing.hpp>
#include <Lielab/utils.hpp>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/src/QR/ColPivHouseholderQR.h>

#include <functional>
#include <limits>

namespace Lielab::optimize
{

NewtonMinimize::NewtonMinimize()
{

}

Eigen::VectorXd NewtonMinimize::operator()(const EuclideanExtremizationSystem system, const Eigen::VectorXd& guess, const ExtremizationOptions options)
{
    this->iteration = 0;
    this->success = false;
    this->message = "";

    ExtremizationOptions _options = options;

    Eigen::VectorXd xcurrent = guess;
    double fcurrent = system.objective(xcurrent);

    if (std::isnan(fcurrent))
    {
        this->message = "NewtonMinimize: NaN found in objective function.";
        this->success = false;
        return xcurrent;
    }

    if (std::isinf(fcurrent))
    {
        this->message = "NewtonMinimize: Inf found in objective function.";
        this->success = false;
        return xcurrent;
    }

    Eigen::VectorXd gradient = system.jacobian(guess);

    const bool gradient_is_defined = !(gradient.size() == 1 && std::isnan(gradient(0, 0)));
    lielab_assert(gradient_is_defined, "Jacobian function must be defined.");

    Eigen::MatrixXd hessian = system.hessian(guess);

    const bool hessian_is_defined = !(hessian.rows() == 1 && hessian.cols() == 2 && std::isnan(hessian(0, 0)) && std::isnan(hessian(0, 1)));
    lielab_assert(hessian_is_defined, "Hessian function must be defined.");

    lielab_assert(gradient.size() == hessian.rows() && gradient.size() == hessian.cols(), "Jacobian must have same size as Hessian. Jacobian was: (" + std::to_string(gradient.size()) + "). Hessian was: (" + std::to_string(hessian.rows()) + ", " + std::to_string(hessian.cols()) + ").");

    const ptrdiff_t n_lower = system.lower_bound.size();
    const ptrdiff_t n_upper = system.upper_bound.size();

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

    // LineSearch search = LineSearch();
    // double alpha = options.initial_alpha;
    // int no_progress = 0;
    // int max_no_progress = 1;
    // double mult_no_progress = 5.0;
    // int no_contractions = 0;
    // int max_no_contractions = 2;

    while (this->success == false)
    {
        // Pre-iteration termination conditions
        if ((_options.max_iterations != -999) && iteration >= _options.max_iterations)
        {
            // End the process if we exceed max iterations. Ignored when max_iterations is -999.
            this->message = "NewtonMinimize: Max iterations reached.";
            this->success = false;
            return xcurrent;
        }

        fcurrent = system.objective(xcurrent);

        if (std::isnan(fcurrent))
        {
            this->message = "NewtonMinimize: NaNs in objective.";
            this->success = false;
            return guess;
        }

        if (std::isinf(fcurrent))
        {
            this->message = "NewtonMinimize: Infs in objective.";
            this->success = false;
            return guess;
        }

        gradient = system.jacobian(xcurrent);

        if (gradient.array().isNaN().any())
        {
            this->message = "NewtonMinimize: NaNs in jacobian.";
            this->success = false;
            return guess;
        }

        if (gradient.array().isInf().any())
        {
            this->message = "NewtonMinimize: Infs in jacobian.";
            this->success = false;
            return guess;
        }

        hessian = system.hessian(xcurrent);

        if (hessian.array().isNaN().any())
        {
            this->message = "NewtonMinimize: NaNs in hessian.";
            this->success = false;
            return guess;
        }

        if (hessian.array().isInf().any())
        {
            this->message = "NewtonMinimize: Infs in hessian.";
            this->success = false;
            return guess;
        }
        
        // Get hessian eigenvalues
        Eigen::EigenSolver<Eigen::MatrixXd> solver(hessian);

        if (solver.info() != Eigen::Success)
        {
            this->message = "NewtonMinimize: Unable to solve for eigenvalues of hessian.";
            this->success = false;
            return guess;
        }

        Eigen::VectorXcd eigenvalues = solver.eigenvalues();

        if ((gradient.cwiseAbs().array() < options.abstol).all())
        {
            // End the method

            bool pos_definite = true;
            bool pos_semidefinite = true;

            for (int ii = 0; ii < eigenvalues.size(); ii++)
            {
                if (eigenvalues(ii).real() < 0.0) pos_definite = false;
                if (eigenvalues(ii).real() <= 0.0 + options.abstol) pos_semidefinite = false;
            }

            if (pos_definite)
            {
                this->message = "NewtonMinimize: Converged to a local minimum.";
                this->success = true;
                return xcurrent;
            }

            if (pos_semidefinite)
            {
                this->message = "NewtonMinimize: Converged, but possibly at a saddle point.";
                this->success = true;
                return xcurrent;
            }

            this->message = "NewtonMinimize: Jacobian too small.";
            this->success = false;
            return xcurrent;
        }

        // Do the Newton iteration.
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(hessian);
        // const int hessian_rank = QR.rank();
        // Eigen::VectorXd xdeltabar;
        // if (hessian_rank == hessian.rows())
        // {
        //     Eigen::PartialPivLU<Eigen::MatrixXd> LU(hessian);
        //     xdeltabar = -LU.solve(gradientbar);
        // }
        // else
        // {
        //     xdeltabar = -QR.solve(gradientbar);
        // }

        const Eigen::VectorXd xdelta = -QR.solve(gradient);

        // TODO: Check LU results and maybe use pseudoinverse? or QR?

        Eigen::VectorXd xnext = xcurrent + xdelta; // TODO: Use line search here instead
        // Use the LineSearch to come up with the next point
        // const double magnitude = system.inner(xcurrent, xdelta, xdelta);

        // Use the alpha the previous iteration has recommended
        // _options.initial_alpha = alpha;
        // CompositeManifold xnext = search(system, xcurrent, xdelta, std::pow(magnitude, 1.0), _options);

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

        const Eigen::VectorXd diff = xnext - xcurrent;
        
        if ((diff.cwiseAbs().array() < _options.abstol*0.1).all())
        {
            // End if progress is not being made.
            this->message = "NewtonMinimize: Progress is not being made (change in x too small).";
            this->success = false;
            return xnext;
        }

        // Check if progress was made
        // const double fnext = system.objective(xnext);
        // if (fnext >= fcurrent)
        // {
        //     // Increase counter of no progress made
        //     no_progress++;
        // }
        // else
        // {
        //     // Progress made, reset counter.
        //     no_progress = 0;
        // }

        // if (no_progress >= max_no_progress)
        // {
        //     this->success = false;
        //     // this->status = ;
        //     this->message = "NewtonMinimize: Objective function is no longer decreasing.";
        //     return xcurrent;
        // }

        // if (fnext < fcurrent)
        // {
        //     // Progress made, reuse the alpha that worked.
        //     alpha = search.alpha;

        //     if (search.iteration == 1)
        //     {
        //         // LineSearch didn't change the search size (no contractions used)
        //         no_contractions++;

        //         if (no_contractions > max_no_contractions)
        //         {
        //             // LineSearch didn't change the search size multiple iterations in a row.
        //             // Increase search size to speed up the optimization.
        //             alpha /= _options.contraction_factor;
        //         }
        //     }
        //     else
        //     {
        //         // LineSearch changed search size. Reset our counter.
        //         no_contractions = 0;
        //     }
        // }
        // else
        // {
        //     // If no progress made, increase the search range.
        //     alpha *= mult_no_progress;
        // }

        // Set up for the next iteration
        this->iteration++;
        xcurrent = xnext;
        // fcurrent = fnext;
    }

    // This shouldnt be called.
    this->message = "NewtonMinimize: Unknown error.";
    this->success = false;
    return guess;
}

}
