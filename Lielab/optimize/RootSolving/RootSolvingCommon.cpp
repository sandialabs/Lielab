#include "RootSolvingCommon.hpp"

#include <functional>
#include <iostream>
namespace Lielab::optimize
{

EuclideanRootSystem::EuclideanRootSystem()
{
    
}

EuclideanRootSystem::EuclideanRootSystem(EuclideanRootSystem_objective_t objective_)
{
    this->objective = objective_;
}

EuclideanRootSystem wrap_with_finite_difference(const EuclideanRootSystem& system, const Eigen::VectorXd& x, const double dx)
{
    const int n_vars = static_cast<int>(x.size());

    EuclideanRootSystem wrapped_sys(system.objective);

    // TODO: Check objective here
    const Eigen::VectorXd _j0 = wrapped_sys.objective(x);
    const int n_funcs = static_cast<int>(_j0.size());

    // Handle the objective jacobian
    const Eigen::MatrixXd djdx = system.jacobian(x);

    if (djdx.rows() == 1 && djdx.cols() == 2 && std::isnan(djdx(0, 0)) && std::isnan(djdx(0, 1)))
    {
        auto wrapped_jacobian = [n_funcs, n_vars, system, dx](const Eigen::VectorXd& _x) -> Eigen::MatrixXd
        {
            Eigen::VectorXd x_copy = _x;
            Eigen::MatrixXd _djdx(n_funcs, n_vars);

            const Eigen::VectorXd j0 = system.objective(_x);

            for (int ii = 0; ii < n_vars; ii++)
            {
                const double x_saved = _x(ii);
                x_copy(ii) += dx;
                const double _dx = x_copy(ii) - x_saved;
                const Eigen::VectorXd j1 = system.objective(x_copy);
                for (int jj = 0; jj < n_funcs; jj++)
                {
                    _djdx(jj, ii) = (j1(jj) - j0(jj))/_dx;
                }
                x_copy(ii) = x_saved;
            }

            return _djdx;
        };

        wrapped_sys.jacobian = wrapped_jacobian;
    }
    else
    {
        wrapped_sys.jacobian = system.jacobian;
    }

    // Handle the objective hessian
    // TODO:

    return wrapped_sys;
}

}
