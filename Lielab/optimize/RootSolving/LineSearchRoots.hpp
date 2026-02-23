#ifndef LIELAB_OPTIMIZE_ROOTSOLVING_LINESEARCHROOTS_HPP
#define LIELAB_OPTIMIZE_ROOTSOLVING_LINESEARCHROOTS_HPP

#include "RootSolvingCommon.hpp"

#include <Eigen/Core>

#include <functional>
#include <string>

namespace Lielab::optimize
{

class LineSearchRoots
{
    public:

    // Metadata
    int num_objective_evals = 0;
    int iteration = 0;
    bool success = false;
    std::string message = "";

    LineSearchRoots();
    Eigen::VectorXd operator()(const EuclideanRootSystem system, const Eigen::VectorXd& guess, const Eigen::VectorXd& direction, const double magnitude, const RootSolvingOptions options);

};

}

#endif
