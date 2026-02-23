#ifndef LIELAB_OPTIMIZE_ROOTSOLVING_HYBRIDROOTSEARCH_HPP
#define LIELAB_OPTIMIZE_ROOTSOLVING_HYBRIDROOTSEARCH_HPP

#include "RootSolvingCommon.hpp"

#include <Eigen/Core>

#include <functional>
#include <string>

namespace Lielab::optimize
{

class HybridRootSearch
{
    public:

    // Metadata
    // int num_objective_evals = 0;
    int iterations = 0;
    bool success = false;
    std::string message = "";

    HybridRootSearch();
    Eigen::VectorXd operator()(const EuclideanRootSystem system, const Eigen::VectorXd& x_guess, const RootSolvingOptions options);

};

}

#endif
