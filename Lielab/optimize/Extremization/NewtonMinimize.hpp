#ifndef LIELAB_OPTIMIZE_ROOTSOLVING_NEWTONMINIMIZE_HPP
#define LIELAB_OPTIMIZE_ROOTSOLVING_NEWTONMINIMIZE_HPP

#include "ExtremizationCommon.hpp"

#include <Eigen/Core>

#include <functional>
#include <string>

namespace Lielab::optimize
{

class NewtonMinimize
{
    public:

    // Metadata
    int iteration = 0;
    bool success = false;
    std::string message = "";

    NewtonMinimize();
    Eigen::VectorXd operator()(const EuclideanExtremizationSystem system, const Eigen::VectorXd& guess, const ExtremizationOptions options);

};

}

#endif
