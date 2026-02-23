#ifndef LIELAB_OPTIMIZE_ROOTSOLVING_NEWTONROOTSEARCH_HPP
#define LIELAB_OPTIMIZE_ROOTSOLVING_NEWTONROOTSEARCH_HPP

#include "RootSolvingCommon.hpp"

#include <Eigen/Core>

#include <functional>
#include <string>

namespace Lielab::optimize
{

class NewtonRootSearch
{
    public:

    // Metadata
    int num_objective_evals = 0;
    int iteration = 0;
    bool success = false;
    std::string message = "";

    NewtonRootSearch();
    Eigen::VectorXd operator()(const EuclideanRootSystem system, const Eigen::VectorXd& x_guess, const RootSolvingOptions options);

};

}

#endif
