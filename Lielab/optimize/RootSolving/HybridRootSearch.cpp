#include "HybridRootSearch.hpp"
#include "RootSolvingCommon.hpp"

#include <Lielab/utils.hpp>

#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>

#include <cassert>
#include <functional>
#include <limits>

#include <iostream>

namespace Lielab::optimize
{

template <typename Field>
struct EuclideanRootFunctor
{
    typedef Eigen::Matrix<Field, Eigen::Dynamic, 1> InputType;
    typedef Eigen::Matrix<Field, Eigen::Dynamic, 1> ValueType;
    typedef Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
    typedef Field Scalar;
    enum { InputsAtCompileTime = Eigen::Dynamic, ValuesAtCompileTime = Eigen::Dynamic };
    int n_variables;
    int n_functions;

    EuclideanRootSystem system;

    EuclideanRootFunctor(const int variables, const int functions)
    {
        this->n_variables = variables;
        this->n_functions = functions;
    }

    int inputs() const
    {
        return this->n_variables;
    }

    int values() const
    {
        return this->n_functions;
    }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& funcs) const
    {
        assert(x.size() == this->n_variables);
        assert(funcs.size() == this->n_functions);

        // TODO: Edit in place for better peformance
        funcs = this->system.objective(x);

        assert(funcs.size() == this->n_functions);

        if (funcs.array().isNaN().any()) return -1;
        if (funcs.array().isInf().any()) return -1;

        return 0;
    }

    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian) const
    {
        assert(x.size() == this->n_variables);
        assert(jacobian.rows() == this->n_variables);
        assert(jacobian.cols() == this->n_functions);

        jacobian = this->system.jacobian(x).transpose();

        assert(jacobian.rows() == this->n_variables);
        assert(jacobian.cols() == this->n_functions);

        if (jacobian.array().isNaN().any()) return -1;
        if (jacobian.array().isInf().any()) return -1;

        return 0;
    }

};

HybridRootSearch::HybridRootSearch()
{

}

Eigen::VectorXd HybridRootSearch::operator()(const EuclideanRootSystem system, const Eigen::VectorXd& x_guess, const RootSolvingOptions options)
{
    this->message = "";
    this->success = false;

    if (system.lower_bound.size() != 0 || system.upper_bound.size() != 0)
    {
        this->message = "HybridRootSearch: lower_bound and upper_bound not implemented.";
        this->success = false;
        return x_guess;
    }

    const Eigen::VectorXd f0 = system.objective(x_guess);

    if (f0.size() != x_guess.size())
    {
        this->message = "HybridRootSearch: Not implemented for under/over-determined (num_functions != num_variables) systems.";
        this->success = false;
        return x_guess;
    }

    const Eigen::MatrixXd dfdx0 = system.jacobian(x_guess);

    if (dfdx0.rows() == 1 && dfdx0.cols() == 2 && std::isnan(dfdx0(0, 0)) && std::isnan(dfdx0(0, 1)))
    {
        this->success = false;
        this->message = "HybridRootSearch: Jacobian function must be defined.";
        return x_guess;
    }

    if (dfdx0.rows() != dfdx0.cols())
    {
        this->success = false;
        this->message = "HybridRootSearch: Jacobian function must be square (got n_rows != n_cols).";
        return x_guess;
    }

    const int n_variables = static_cast<int>(x_guess.size());
    const int n_functions = static_cast<int>(f0.size());

    Eigen::VectorXd x = x_guess;

    // Set up the optimization routine
    EuclideanRootFunctor<double> functor(n_variables, n_functions);
    functor.system = system;
    // Eigen::NumericalDiff<EuclideanRootFunctor<double>> _functor(functor);
    // Eigen::HybridNonLinearSolver<Eigen::NumericalDiff<EuclideanRootFunctor<double>>> solver(_functor);

    Eigen::HybridNonLinearSolver<EuclideanRootFunctor<double>> solver(functor);
    solver.parameters.xtol = options.tol;

    // solver.diag.setConstant(n, 1.0); // Scaling parameters. Probably not needed, just let the user scale on their own.
    // solver.useExternalScaling = true;

    // Do the optimization
    const int info = solver.solve(x);
    
    // Postprocess
    this->iterations = static_cast<int>(solver.iter);

    // Info could be:
    // Running = -1
    // ImproperInputParameters = 0
    if (info == 1)
    {
        // RelativeErrorTooSmall = 1
        this->success = true;
        this->message = "HybridRootSearch: Converged to tolerance.";
    }
    else if (info == 2)
    {
        // TooManyFunctionEvaluation = 2
        this->success = false;
        this->message = "HybridRootSearch: Too many function evaluations.";
    }
    else if (info == 3)
    {
        // TolTooSmall = 3
        this->success = false;
        this->message = "HybridRootSearch: Tolerance too small.";
    }
    else if (info == 4)
    {
        // NotMakingProgressJacobian = 4
        this->success = false;
        this->message = "HybridRootSearch: Not making progress.";
    }
    else if (info == 5)
    {
        // NotMakingProgressIterations = 5
        this->success = false;
        this->message = "HybridRootSearch: x is unchanging.";
    }
    else if (info == 6)
    {
        // UserAsked = 6
        this->success = false;
        this->message = "HybridRootSearch: NaNs or Infs in objective or jacobian.";
    }

    return x;
}

}
