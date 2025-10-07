#ifndef LIELAB_INTEGRATE_ODESOLUTION_HPP
#define LIELAB_INTEGRATE_ODESOLUTION_HPP

#include "Lielab/domain.hpp"

#include <Eigen/Core>

#include <limits>
#include <vector>

namespace Lielab::integrate
{

/*!

*/
class ODESolution
{
    public:

    // Solution metadata
    bool success = false;
    int status = std::numeric_limits<int>::quiet_NaN();
    std::string message = "";
    double time_to_solution = std::numeric_limits<double>::quiet_NaN();
    Eigen::VectorXd residuals = Eigen::VectorXd::Zero(0);
    Eigen::MatrixXd debug;

    // ODE path data
    Eigen::VectorXd t = Eigen::VectorXd::Zero(0);
    std::vector<Lielab::domain::CompositeManifold> y;
    Eigen::MatrixXd ybar = Eigen::MatrixXd::Zero(0, 0);
    std::vector<Lielab::domain::CompositeAlgebra> theta;
    Eigen::MatrixXd thetabar = Eigen::MatrixXd::Zero(0, 0);
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(0);
    Eigen::VectorXd p = Eigen::VectorXd::Zero(0);
    
    // Chunking algorithm params
    size_t chunk_size = 2048;
    size_t current_index = 0;

    ODESolution();
    ODESolution(const ODESolution& other);
    ODESolution(const size_t num_eoms);
    ODESolution& operator=(const ODESolution& other);
    ODESolution copy() const;
    void trim_chunk(const size_t last_index);
    void add_chunk();
    void add_data(const double t_add, const Eigen::VectorXd& ybar_add);
    void add_data(const double t_add, const Eigen::VectorXd& ybar_add, const Eigen::VectorXd& thetabar_add);
};

}

#endif
