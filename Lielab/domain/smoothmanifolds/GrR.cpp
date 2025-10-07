#include "GrR.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

namespace Lielab::domain
{

std::string GrR::to_string() const
{
    const std::string k_str = std::to_string(this->data.cols());
    const std::string n_str = std::to_string(this->data.rows());

    return "Gr(" + k_str + ", " + n_str + ", R)";
}

GrR::GrR()
{
    this->data = Eigen::MatrixXd::Identity(0, 0);
}

GrR::GrR(const size_t k, const size_t n)
{
    if (k > n)
    {
        throw Lielab::utils::Error("n must be greater than or equal to k.");
    }

    this->data = Eigen::MatrixXd::Identity(n, k);
}

size_t GrR::get_dimension() const
{
    const size_t k = this->data.cols();
    const size_t n = this->data.rows();
    return k*(n - k);
}

size_t GrR::get_size() const
{
    const size_t k = this->data.cols();
    const size_t n = this->data.rows();
    return k*n;
}

double GrR::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
{
    /*!
    *
    * Gets a value in the matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t k = this->data.cols();
    const size_t n = this->data.rows();

    if (n == 0) return nan;
    if (k == 0) return nan;

    if (index1 >= static_cast<ptrdiff_t>(n)) return nan;
    if (index2 >= static_cast<ptrdiff_t>(k)) return nan;
    if (std::abs(index1) > static_cast<ptrdiff_t>(n)) return nan;
    if (std::abs(index2) > static_cast<ptrdiff_t>(k)) return nan;

    size_t _index1;
    if (index1 < 0)
    {
        _index1 = static_cast<size_t>(static_cast<ptrdiff_t>(n) + index1);
    }
    else
    {
        _index1 = static_cast<size_t>(index1);
    }

    size_t _index2;
    if (index2 < 0)
    {
        _index2 = static_cast<size_t>(static_cast<ptrdiff_t>(k) + index2);
    }
    else
    {
        _index2 = static_cast<size_t>(index2);
    }

    return this->data(_index1, _index2);
}

Eigen::VectorXd GrR::serialize() const
{
    /*!
    * 
    * Returns a serialized representation.
    */

    return this->data.reshaped<Eigen::RowMajor>();
}

void GrR::unserialize(const Eigen::VectorXd& vec)
{
    /*!
    * 
    * Sets the GrR object from a serialized vector.
    */

    const size_t vdim = vec.size();
    const size_t cols = this->data.cols();
    const size_t rows = this->data.rows();

    const size_t max_ind = std::min(cols*rows, vdim);

    for (size_t vind = 0; vind < max_ind; vind++)
    {
        const size_t row = static_cast<size_t>(std::floor(vind / cols));
        const size_t col = vind % cols;
        this->data(row, col) = vec(vind);
    }
}

void GrR::unserialize(std::initializer_list<double> vec)
{
    this->unserialize(Eigen::VectorXd{std::move(vec)});
}

Eigen::MatrixXd GrR::get_matrix() const
{
    return this->data;
}

Eigen::VectorXd GrR::project_onto(const Eigen::VectorXd& vec) const
{
    // TODO: Linguistically, this function name is backwards but project() is already taken.
    // TODO: Some of these quantities can probably be computed at object initialization.

    return this->data*(this->data.transpose()*this->data).inverse()*this->data.transpose()*vec;
}

}
