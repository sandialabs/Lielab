#include "SP.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

SP::SP() : SP(0)
{
    /*! \f{equation*}{ () \rightarrow SP \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::SP x, y, z;
    * 
    */

}

// SP::~SP();

SP::SP(const SP::matrix_t& matrix)
{
    /*! \f{equation}{(\mathbb{R}^{n \times n}) \rightarrow SP \f}
    *
    * Constructor instantiating an \f$SP\f$ object from an
    * \f$n \times n\f$ real matrix.
    *
    */

    lielab_assert(matrix.rows() % 2 == 0, "Input matrix shape must be even.");
    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->point.noalias() = matrix;
}

SP SP::identity(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SP \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The SP element. 
    */

    return SP(shape);
}

// TODO: Projection

SP::SP(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SP \f}
    *
    * Constructor instantiating an \f$SP\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::SP x(2), y(4), z(6);
    * 
    * @param[in] shape The shape of the data matrix.
    */

    lielab_assert(shape % 2 == 0, "Shape of SP must be even.");

    this->point.noalias() = Eigen::MatrixXd::Identity(shape, shape);
}

std::string SP::to_string() const
{
    return "SP(" + std::to_string(this->get_shape()) + ", R)";
}

int SP::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    return this->get_shape() * (this->get_shape() + 1) / 2;
}

int SP::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    return static_cast<int>(std::pow(this->get_shape(), 2));
}

bool SP::is_abelian() const
{
    return false;
}

int SP::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return static_cast<int>(this->point.rows());
}

SP::point_t SP::get_point() const
{
    return this->point;
}

Eigen::VectorXd SP::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */

    return this->point.reshaped<Eigen::RowMajor>();
}

void SP::unserialize(const Eigen::VectorXd& serialized)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the SP object from a serialized vector.
    */

    const int vdim = static_cast<int>(serialized.size());
    const int max_ind = std::min(this->get_size(), vdim);

    for (int vind = 0; vind < max_ind; vind++)
    {
        const int row = vind / this->get_shape();
        const int col = vind % this->get_shape();
        this->point(row, col) = serialized(vind);
    }
}

void SP::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

SP::matrix_t SP::get_matrix() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times n} \f}
    * 
    * Returns a matrix representation.
    * 
    * Formerly called "get_ados_representation()".
    * 
    * Ado, Igor D. "Note on the representation of finite continuous groups by
    *               means of linear substitutions, Izv. Fiz." Mat. Obsch.(Kazan)
    *               7.1 (1935): 935.
    * 
    * Ado, Igor D. "The representation of Lie algebras by matrices." Uspekhi
    *               Matematicheskikh Nauk 2.6 (1947): 159-173.
    */

    return this->point;
}

SP::field_t SP::operator()(const int index1, const int index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */
    
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int shape = this->get_shape();
    if (shape == 0) return nan;

    // If input index is negative, index from the back of the array
    const int _index1 = (index1 < 0) ? shape + index1 : index1;
    const int _index2 = (index2 < 0) ? shape + index2 : index2;

    // Error check for out of bounds
    if (_index1 < 0) return nan;
    if (_index1 >= shape) return nan;
    if (_index2 < 0) return nan;
    if (_index2 >= shape) return nan;

    return this->point(_index1, _index2);
}

SP SP::operator*(const SP& other) const
{
    /*! \f{equation*}{ (SP, SP) \rightarrow SP \f}
    *
    * Group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return SP(this->point * other.point);
}

SP& SP::operator*=(const SP& other)
{
    /*! \f{equation*}{ (SP, SP) \rightarrow SP \f}
    *
    * In place group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point *= other.point;
    return *this;
}

SP SP::inverse() const
{
    /*! \f{equation*}{ (SP) \rightarrow SP \f}
    * 
    * Returns the inverse.
    */

    return SP(this->point.inverse());
}

}
