#ifndef LIELAB_OPTIMIZE_EXTREMIZATION_GOLDENMINIMIZE_HPP
#define LIELAB_OPTIMIZE_EXTREMIZATION_GOLDENMINIMIZE_HPP

#include "ExtremizationCommon.hpp"

#include <Eigen/Core>

#include <functional>
#include <limits>
#include <string>

namespace Lielab::optimize
{

class GoldenMinimize
{
    public:

    // general
    int iteration = 0;
    bool success = false;
    std::string message = "";

    GoldenMinimize();
    Eigen::VectorXd operator()(EuclideanExtremizationSystem system, const ExtremizationOptions options);
};

}

#endif
