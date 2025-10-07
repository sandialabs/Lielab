#include "SE.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string SE::to_string() const
{
    if (this->_shape == 0) return "SE(0)";
    return "SE(" + std::to_string(this->_shape-1) + ")";
}

SE::SE() : SE(0)
{
    /*! \f{equation*}{ () \rightarrow SE \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::SE x, y, z;
    * 
    */

}

SE::SE(const size_t n)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SE \f}
    *
    * Constructor instantiating an \f$SE\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::SE x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */

    this->_shape = n + 1;
    this->data = Eigen::MatrixXd::Identity(n + 1, n + 1);
}

SE SE::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{SE} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The SE element. 
    */

    if (shape == 0)
    {
        SE out;
        out._shape = 0;
        out.data = Eigen::MatrixXd::Zero(0, 0);
        return out;
    }

    return SE(shape - 1);
}

// Eigen::MatrixXd SE::project(const Eigen::MatrixXd & other)
// {
//     /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{SE} \f}
//     *
//     * Projects a matrix suitable for data.
//     */

//     if (other.rows() != other.cols())
//     {
//         throw Lielab::utils::Error("Size of the matrix must be square.");
//     }

//     Eigen::HouseholderQR<Eigen::MatrixXd> qr(other);
//     Eigen::MatrixXd Q = qr.householderQ();
//     Eigen::MatrixXd R = qr.matrixQR();
//     Eigen::MatrixXd P = Eigen::MatrixXd::Zero(other.rows(), other.cols());

//     for (int ii = 0; ii < other.rows(); ii++)
//     {
//         if (R(ii,ii) < 0)
//         {
//             P(ii, ii) = -1.0;
//         }
//         else
//         {
//             P(ii, ii) = 1.0;
//         }
//     }

//     return Q*P;
// }

size_t SE::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    return (this->_shape - 1) * (this->_shape - 2) / 2 + this->_shape - 1;
}

size_t SE::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return this->_shape;
}

size_t SE::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    return static_cast<size_t>(std::pow(this->_shape, 2));
}

Eigen::MatrixXd SE::get_matrix() const
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

    return this->data;
}

SE SE::inverse() const
{
    /*! \f{equation*}{ (SE) \rightarrow SE \f}
    * 
    * Returns the inverse.
    */

    return this->data.inverse();
}

Eigen::VectorXd SE::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */

    return this->data.reshaped<Eigen::RowMajor>();
}

void SE::unserialize(const Eigen::VectorXd& vec)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the SE object from a serialized vector.
    */

    const size_t vdim = vec.size();
    const size_t max_ind = std::min(this->get_size(), vdim);

    for (size_t vind = 0; vind < max_ind; vind++)
    {
        const size_t row = static_cast<size_t>(std::floor(vind / this->_shape));
        const size_t col = vind % this->_shape;
        this->data(row, col) = vec(vind);
    }
}

void SE::unserialize(std::initializer_list<double> vec)
{
    this->unserialize(Eigen::VectorXd{std::move(vec)});
}

double SE::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t shape = this->get_shape();
    if (shape == 0) return nan;

    if (index1 >= static_cast<ptrdiff_t>(shape)) return nan;
    if (index2 >= static_cast<ptrdiff_t>(shape)) return nan;
    if (std::abs(index1) > static_cast<ptrdiff_t>(shape)) return nan;
    if (std::abs(index2) > static_cast<ptrdiff_t>(shape)) return nan;

    size_t _index1;
    if (index1 < 0)
    {
        _index1 = static_cast<size_t>(static_cast<ptrdiff_t>(shape) + index1);
    }
    else
    {
        _index1 = static_cast<size_t>(index1);
    }

    size_t _index2;
    if (index2 < 0)
    {
        _index2 = static_cast<size_t>(static_cast<ptrdiff_t>(shape) + index2);
    }
    else
    {
        _index2 = static_cast<size_t>(index2);
    }

    return this->data(_index1, _index2);
}

SE SE::operator*(const SE & other) const
{
    /*! \f{equation*}{ (SE, SE) \rightarrow SE \f}
    *
    * Group product.
    */

    assert(this->_shape == other.get_shape());
    return SE(this->data * other.data);
}

SE & SE::operator*=(const SE & other)
{
    /*! \f{equation*}{ (SE, SE) \rightarrow SE \f}
    *
    * In place group product.
    */

    assert(this->_shape == other.get_shape());
    this->data *= other.data;
    return *this;
}

std::ostream & operator<<(std::ostream& os, const SE & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << static_cast<const Eigen::MatrixXd>(other.data);
    return os;
}

}
