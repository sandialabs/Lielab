#ifndef LIELAB_OPTIMIZE_EXTREMIZATION_LINESEARCH_HPP
#define LIELAB_OPTIMIZE_EXTREMIZATION_LINESEARCH_HPP

#include "ExtremizationCommon.hpp"

#include <Eigen/Core>

#include <functional>
#include <limits>
#include <string>

namespace Lielab::optimize
{

class LineSearch
{
    public:

    // General
    int iteration = 0;
    bool success = false;
    std::string message = "";

    // LineSearch specific
    double alpha = std::numeric_limits<double>::quiet_NaN();

    LineSearch();
    Eigen::VectorXd operator()(const EuclideanExtremizationSystem system,
                               const Eigen::VectorXd& guess,
                               const Eigen::VectorXd& direction,
                               const double magnitude,
                               const ExtremizationOptions options);
};

}

#endif
