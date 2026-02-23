#include "Grassmannian.hpp"

#include "Lielab/domain/liealgebras/rn.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/QR/ColPivHouseholderQR.h>

namespace Lielab::domain
{

Grassmannian::Grassmannian()
{
    this->point.noalias() = Eigen::VectorXd::Zero(0);
    this->axes.noalias() = Eigen::MatrixXd::Identity(0, 0);
    this->point_projector.noalias() = Eigen::MatrixXd::Zero(0, 0);
}

Grassmannian::Grassmannian(const int k, const int n)
{
    lielab_assert(k <= n, "n (" + std::to_string(n) + ") must be greater than or equal to k (" + std::to_string(k) + ").");

    this->point.noalias() = Eigen::VectorXd::Zero(n);
    this->axes.noalias() = Eigen::MatrixXd::Identity(n, k);
    this->point_projector.noalias() = this->axes*(this->axes.transpose()*this->axes).inverse()*this->axes.transpose();
}

Grassmannian::Grassmannian(const Eigen::VectorXd& point, const Eigen::MatrixXd& axes)
{
    lielab_assert(axes.cols() <= axes.rows(), "Axes rows (" + std::to_string(axes.rows()) + ") must be greater than or equal to axes cols (" + std::to_string(axes.cols()) + ").");
    lielab_assert(point.size() == axes.rows(), "Size of point (" + std::to_string(point.size()) + ") must be the same as axes rows (" + std::to_string(axes.rows()) + ").");

    this->point.noalias() = point;
    this->axes.noalias() = axes;
    this->point_projector.noalias() = this->axes*(this->axes.transpose()*this->axes).inverse()*this->axes.transpose();
}

Grassmannian Grassmannian::project(const Eigen::VectorXd& point, const Eigen::MatrixXd& axes)
{
    lielab_assert(axes.cols() <= axes.rows(), "Axes rows (" + std::to_string(axes.rows()) + ") must be greater than or equal to axes cols (" + std::to_string(axes.cols()) + ").");
    lielab_assert(point.size() == axes.rows(), "Size of point (" + std::to_string(point.size()) + ") must be the same as axes rows (" + std::to_string(axes.rows()) + ").");

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(axes);
    const int rank = static_cast<int>(QR.rank());

    lielab_assert(rank == axes.cols(), "Rank of axis must equal number of columns.");

    const Eigen::MatrixXd Q = QR.householderQ();
    const Eigen::MatrixXd orth_axes = Q.block(0, 0, axes.rows(), axes.cols());

    Grassmannian out(point, orth_axes);
    out.point.noalias() = out.project_point(point);
    return out;
}

std::string Grassmannian::to_string() const
{
    const std::string k_str = std::to_string(this->axes.cols());
    const std::string n_str = std::to_string(this->axes.rows());

    return "Grassmannian(" + k_str + ", " + n_str + ", R)";
}

int Grassmannian::get_dimension() const
{
    return static_cast<int>(this->axes.cols()); // k
}

int Grassmannian::get_size() const
{
    return static_cast<int>(this->axes.rows()); // n
}

Grassmannian::point_t Grassmannian::get_point() const
{
    return this->point;
}

Eigen::VectorXd Grassmannian::serialize() const
{
    /*!
    * 
    * Returns a serialized representation.
    */

    return this->point.reshaped<Eigen::RowMajor>();
}

void Grassmannian::unserialize(const Eigen::VectorXd& vec)
{
    /*!
    * 
    * Sets the Grassmannian object from a serialized vector.
    */

    const int sz = std::min(this->get_size(), static_cast<int>(vec.size()));
    this->point(Eigen::seqN(0, sz)) = vec(Eigen::seqN(0, sz));
}

void Grassmannian::unserialize(std::initializer_list<double> vec)
{
    this->unserialize(Eigen::VectorXd{std::move(vec)});
}

const double& Grassmannian::operator[](const int index) const
{
    const int len = static_cast<int>(this->point.size());

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? len + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for Grassmannian of size " + std::to_string(len));

    return this->point(_index);
}

double& Grassmannian::operator[](const int index)
{
    const int len = static_cast<int>(this->point.size());

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? len + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for Grassmannian of size " + std::to_string(len));

    return this->point(_index);
}

Eigen::VectorXd Grassmannian::project_point(const Eigen::VectorXd& other_point) const
{
    lielab_assert(other_point.size() == this->axes.rows(), "Size of point (" + std::to_string(other_point.size()) + ") must be the same as axes rows (" + std::to_string(this->axes.rows()) + ").");

    return this->point_projector*other_point;
}


Eigen::VectorXd Grassmannian::project_vector_onto_tangent_space(const Eigen::VectorXd& vector) const
{
    lielab_assert(vector.size() == this->axes.rows(), "Size of vector (" + std::to_string(vector.size()) + ") must be the same as the total space (" + std::to_string(this->axes.rows()) + ").");
    return this->axes*(this->axes.transpose()*vector);
}

Eigen::VectorXd Grassmannian::project_vector_onto_normal_space(const Eigen::VectorXd& vector) const
{
    lielab_assert(vector.size() == this->axes.rows(), "Size of vector (" + std::to_string(vector.size()) + ") must be the same as the total space (" + std::to_string(this->axes.rows()) + ").");
    return vector - this->project_vector_onto_tangent_space(vector);
}

Eigen::VectorXd Grassmannian::axes_intersection(const Eigen::VectorXd& other_point, const Eigen::MatrixXd& other_axes) const
{
    using Lielab::utils::horizontal_stack;

    lielab_assert(other_point.size() == this->axes.rows(), "Other point must have same size as total space dimension.");
    lielab_assert(this->axes.rows() == other_axes.rows(), "Both axes must have same number of rows.");
    lielab_assert(this->axes.cols() + other_axes.cols() == this->axes.rows(), "Combined axes must equal n-space."); // TODO: Only point-wise intersections for now (no Grassmannian subspaces)

    const Eigen::MatrixXd combined_axes = horizontal_stack<double>({this->axes, -other_axes});

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(combined_axes);
    const int rank = static_cast<int>(QR.rank());

    // Rank-deficient. Doesn't intersect. Return nans.
    if (rank < this->axes.rows()) return std::numeric_limits<double>::quiet_NaN()*other_point;

    const Eigen::VectorXd sol = QR.solve(other_point);

    return this->axes*sol(Eigen::seqN(0, this->axes.cols()));
}

}
