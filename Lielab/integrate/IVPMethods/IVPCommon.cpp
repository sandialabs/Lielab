#include "IVPCommon.hpp"

#include <string>

namespace Lielab::integrate
{

std::string IVPSolution::to_string() const
{
    std::string out;
    const std::string success_str = (this->success) ? "True" : "False";

    out += "  success: " + success_str + "\n";
    out += "   status: " + std::to_string(static_cast<int>(this->status)) + ": " + this->message + "\n";
    out += "        t: Array (" + std::to_string(this->t.size()) + ",)\n";
    out += "        y: CompositeManifold (" + std::to_string(this->y.size()) + ",)\n";
    return out;
}

IVPSolution::IVPSolution()
{
    // this->t = Eigen::VectorXd::Zero(this->chunk);
    // this->ybar = Eigen::MatrixXd::Zero(this->chunk, 0);
    // this->thetabar = Eigen::MatrixXd::Zero(this->chunk, 0);
}

IVPSolution::IVPSolution(const IVPSolution& other)
{
    /*!
    * Copy constructor for IntegralCurve
    */

    this->message = other.message;
    this->status = other.status;
    this->success = other.success;
    this->time_to_solution = other.time_to_solution;
    
    this->chunk_size = other.chunk_size;
    this->current_index = other.current_index;

    this->t = other.t;
    this->y = other.y;
    this->ybar = other.ybar;
    this->theta = other.theta;
    this->thetabar = other.thetabar;

    this->debug = other.debug;
}

IVPSolution::IVPSolution(const size_t num_eoms)
{
    this->t = Eigen::VectorXd::Zero(0);
    this->ybar = Eigen::MatrixXd::Zero(0, num_eoms);
    this->thetabar = Eigen::MatrixXd::Zero(0, num_eoms);
}

IVPSolution& IVPSolution::operator=(const IVPSolution& other)
{
    this->message = other.message;
    this->status = other.status;
    this->success = other.success;
    this->time_to_solution = other.time_to_solution;
    
    this->chunk_size = other.chunk_size;
    this->current_index = other.current_index;

    this->t = other.t;
    this->y = other.y;
    this->ybar = other.ybar;
    this->theta = other.theta;
    this->thetabar = other.thetabar;

    this->debug = other.debug;
    return *this;
}

void IVPSolution::trim_chunk()
{
    const size_t sz_manifold = this->ybar.cols();
    const size_t sz_algebra = this->thetabar.cols();

    this->t.conservativeResize(this->current_index);
    this->ybar.conservativeResize(this->current_index, sz_manifold);
    this->thetabar.conservativeResize(this->current_index, sz_algebra);
}

void IVPSolution::trim_chunk(const ptrdiff_t last_index)
{
    // TODO: Is this method needed anymore?
    const ptrdiff_t sz_manifold = this->ybar.cols();
    const ptrdiff_t sz_algebra = this->thetabar.cols();

    this->t.conservativeResize(last_index + 1);
    this->ybar.conservativeResize(last_index + 1, sz_manifold);
    this->thetabar.conservativeResize(last_index + 1, sz_algebra);
}

void IVPSolution::add_chunk()
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

void IVPSolution::add_data(const double t_add, const Eigen::VectorXd& ybar_add)
{
    const int length = static_cast<int>(this->t.size());

    // Resize storage if necessary
    if (this->current_index >= length)
    {
        this->add_chunk();
    }

    // Save data
    const size_t sz_manifold = this->ybar.cols();

    this->t(this->current_index) = t_add;
    this->ybar.block(this->current_index, 0, 1, sz_manifold).noalias() = ybar_add.transpose();
    
    this->current_index += 1;
}

void IVPSolution::add_data(const double t_add, const Eigen::VectorXd& ybar_add, const Eigen::VectorXd& thetabar_add)
{
    const int length = static_cast<int>(this->t.size());

    // Resize storage if necessary
    if (this->current_index >= length)
    {
        this->add_chunk();
    }

    // Save data
    const size_t sz_manifold = this->ybar.cols();
    const size_t sz_algebra = this->thetabar.cols();

    this->t(this->current_index) = t_add;
    this->ybar.block(this->current_index, 0, 1, sz_manifold).noalias() = ybar_add.transpose();
    this->thetabar.block(this->current_index, 0, 1, sz_algebra).noalias() = thetabar_add.transpose();
    
    this->current_index += 1;
}

EuclideanIVPSystem::EuclideanIVPSystem(EuclideanIVP_vectorfield_t vf)
{
    this->vectorfield = vf;
}

HomogeneousIVPSystem::HomogeneousIVPSystem(HomogeneousIVP_generator_t vf)
{
    this->generator = vf;
}

EuclideanClassicHamiltonianIVPSystem::EuclideanClassicHamiltonianIVPSystem(EuclideanClassicHamiltonianIVP_dVdq_t dVdq_, EuclideanClassicHamiltonianIVP_dTdp_t dTdp_)
{
    this->dVdq = dVdq_;
    this->dTdp = dTdp_;
}

}
