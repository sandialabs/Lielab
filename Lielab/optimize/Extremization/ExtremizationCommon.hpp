#ifndef LIELAB_OPTIMIZE_EXTREMIZATION_EXTREMIZATIONCOMMON_HPP
#define LIELAB_OPTIMIZE_EXTREMIZATION_EXTREMIZATIONCOMMON_HPP

#include "Lielab/domain.hpp"

#include <Eigen/Core>

#include <cmath>
#include <functional>

namespace Lielab::optimize
{

enum class ExtremizationMethod
{
    Undefined = 0,
    LineSearch = 1,
    Newton = 2,
    Golden = 3,
};

struct ExtremizationOptions
{
    public:

    // Extremization meta options
    ExtremizationMethod method = ExtremizationMethod::Undefined;

    // Multi method options
    double abstol = 1e-6;
    double reltol = 1e-6;
    int max_iterations = -999;
    double contraction_factor = 0.5;

    // Golden specific
    double golden_ratio = (std::sqrt(5.0) + 1.0)/2.0;

    // LineSearch specific
    double initial_alpha = std::numeric_limits<double>::quiet_NaN();
    double initial_step_size = 1.0;
    double sufficient_decrease = 0.5;
};

using EuclideanExtremizationSystem_objective_t = std::function<double(const Eigen::VectorXd&)>;
using EuclideanExtremizationSystem_jacobian_t = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;
using EuclideanExtremizationSystem_hessian_t = std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>;

struct EuclideanExtremizationSystem
{
    public:
    Eigen::VectorXd lower_bound = Eigen::VectorXd::Zero(0);
    Eigen::VectorXd upper_bound = Eigen::VectorXd::Zero(0);

    EuclideanExtremizationSystem_objective_t objective = [](const Eigen::VectorXd& x){return std::numeric_limits<double>::quiet_NaN();};
    EuclideanExtremizationSystem_jacobian_t jacobian = [](const Eigen::VectorXd& x){return Eigen::VectorXd::Ones(1)*std::numeric_limits<double>::quiet_NaN();};
    EuclideanExtremizationSystem_hessian_t hessian = [](const Eigen::VectorXd& x){return Eigen::MatrixXd::Ones(1, 2)*std::numeric_limits<double>::quiet_NaN();};

    EuclideanExtremizationSystem(EuclideanExtremizationSystem_objective_t objective_);
};

}

#endif
