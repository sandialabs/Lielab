#ifndef LIELAB_OPTIMIZE_ROOTSOLVINGCOMMON_HPP
#define LIELAB_OPTIMIZE_ROOTSOLVINGCOMMON_HPP

#include <Eigen/Core>

#include <functional>

namespace Lielab::optimize
{

enum class RootSolvingMethod
{
    Undefined = 0,
    Newton = 1,
    LineSearch = 2,
    Hybrid = 3,
};

struct RootSolvingOptions
{
    public:

    // Root solving meta options
    RootSolvingMethod method = RootSolvingMethod::Undefined;

    // Multi method options
    double tol = 1e-8;
    double dx = 1e-8;
    int max_iterations = -999;

    // LineSearch specific
    double contraction_factor = 0.5;
    double initial_alpha = std::numeric_limits<double>::quiet_NaN();
    double initial_step_size = 1.0;
    double sufficient_decrease = 0.5;
};

using EuclideanRootSystem_objective_t = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;
using EuclideanRootSystem_jacobian_t = std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>;

struct EuclideanRootSystem
{
    public:

    EuclideanRootSystem_objective_t objective = [](const Eigen::VectorXd& x){return Eigen::VectorXd::Ones(1)*std::numeric_limits<double>::signaling_NaN();};
    EuclideanRootSystem_jacobian_t jacobian = [](const Eigen::VectorXd& x){return Eigen::MatrixXd::Ones(1, 2)*std::numeric_limits<double>::signaling_NaN();};

    Eigen::VectorXd lower_bound = Eigen::VectorXd::Zero(0);
    Eigen::VectorXd upper_bound = Eigen::VectorXd::Zero(0);

    EuclideanRootSystem();
    EuclideanRootSystem(EuclideanRootSystem_objective_t objective_);
};

EuclideanRootSystem wrap_with_finite_difference(const EuclideanRootSystem& system, const Eigen::VectorXd& x, const double dx = 1e-6);

}

#endif
