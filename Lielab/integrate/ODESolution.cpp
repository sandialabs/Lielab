#include "ODESolution.hpp"

#include <Eigen/Core>

#include <fstream>
#include <iomanip>

namespace Lielab::integrate
{

ODESolution::ODESolution()
{
    // this->t = Eigen::VectorXd::Zero(this->chunk);
    // this->ybar = Eigen::MatrixXd::Zero(this->chunk, 0);
    // this->thetabar = Eigen::MatrixXd::Zero(this->chunk, 0);
}

ODESolution::ODESolution(const ODESolution& other)
{
    /*!
    * Copy constructor for IntegralCurve
    */

    this->message = other.message;
    this->status = other.status;
    this->success = other.success;
    
    this->chunk_size = other.chunk_size;

    this->t = other.t;
    this->y = other.y;
    this->ybar = other.ybar;
    this->theta = other.theta;
    this->thetabar = other.thetabar;
    this->x0 = other.x0;
    this->p = other.p;

    this->residuals = other.residuals;
    this->time_to_solution = other.time_to_solution;

    this->debug = other.debug;
}

ODESolution::ODESolution(const size_t num_eoms)
{
    this->t = Eigen::VectorXd::Zero(0);
    this->ybar = Eigen::MatrixXd::Zero(0, num_eoms);
    this->thetabar = Eigen::MatrixXd::Zero(0, num_eoms);
}

ODESolution& ODESolution::operator=(const ODESolution& other)
{
    this->message = other.message;
    this->status = other.status;
    this->success = other.success;
    
    this->chunk_size = other.chunk_size;

    this->t = other.t;
    this->y = other.y;
    this->ybar = other.ybar;
    this->theta = other.theta;
    this->thetabar = other.thetabar;
    this->x0 = other.x0;
    this->p = other.p;

    this->residuals = other.residuals;
    this->time_to_solution = other.time_to_solution;

    this->debug = other.debug;
    return *this;
}

ODESolution ODESolution::copy() const
{
    /*!
    Awkward method to force a copy of the current object. Mainly for Python.
    */
    const ODESolution other = *this;
    return other;
}

void ODESolution::trim_chunk(const size_t last_index)
{
    const size_t sz_manifold = this->ybar.cols();
    const size_t sz_algebra = this->thetabar.cols();

    Eigen::VectorXd temp_t(last_index + 1);
    Eigen::MatrixXd temp_ybar(last_index + 1, sz_manifold);
    Eigen::MatrixXd temp_thetabar(last_index + 1, sz_algebra);

    temp_t = this->t.head(last_index + 1);
    temp_ybar = this->ybar.block(0, 0, last_index + 1, sz_manifold);
    temp_thetabar = this->thetabar.block(0, 0, last_index + 1, sz_algebra);

    this->t = temp_t;
    this->ybar = temp_ybar;
    this->thetabar = temp_thetabar;
}

void ODESolution::add_chunk()
{
    const size_t sz_manifold = this->ybar.cols();
    const size_t sz_algebra = this->thetabar.cols();
    const size_t length = this->t.size();

    Eigen::VectorXd temp_t = Eigen::VectorXd::Zero(length + this->chunk_size);
    Eigen::MatrixXd temp_ybar = Eigen::MatrixXd::Zero(length + this->chunk_size, sz_manifold);
    Eigen::MatrixXd temp_thetabar = Eigen::MatrixXd::Zero(length + this->chunk_size, sz_algebra);

    temp_t.head(length) = this->t;
    temp_ybar.block(0, 0, length, sz_manifold) = this->ybar;
    temp_thetabar.block(0, 0, length, sz_algebra) = this->thetabar;

    this->t = temp_t;
    this->ybar = temp_ybar;
    this->thetabar = temp_thetabar;
}

void ODESolution::add_data(const double t_add, const Eigen::VectorXd& ybar_add)
{
    const size_t length = this->t.size();

    // Resize storage if necessary
    if (this->current_index >= length)
    {
        this->add_chunk();
    }

    // Save data
    const size_t sz_manifold = this->ybar.cols();

    this->t(this->current_index) = t_add;
    this->ybar.block(this->current_index, 0, 1, sz_manifold) = ybar_add.transpose();
    
    this->current_index += 1;
}

void ODESolution::add_data(const double t_add, const Eigen::VectorXd& ybar_add, const Eigen::VectorXd& thetabar_add)
{
    const size_t length = this->t.size();

    // Resize storage if necessary
    if (this->current_index >= length)
    {
        this->add_chunk();
    }

    // Save data
    const size_t sz_manifold = this->ybar.cols();
    const size_t sz_algebra = this->thetabar.cols();

    this->t(this->current_index) = t_add;
    this->ybar.block(this->current_index, 0, 1, sz_manifold) = ybar_add.transpose();
    this->thetabar.block(this->current_index, 0, 1, sz_algebra) = thetabar_add.transpose();
    
    this->current_index += 1;
}

}
