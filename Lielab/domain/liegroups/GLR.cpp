#include "GLR.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string GLR::to_string() const
{
    return "GL(" + std::to_string(this->_shape) + ", R)";
}

GLR::GLR() : GLR(0)
{
    /*! \f{equation*}{ () \rightarrow GLR \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::GLR x, y, z;
    * 
    */

}

GLR::GLR(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow GLR \f}
    *
    * Constructor instantiating an \f$GLR\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::GLR x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */
    
    this->data = data_t::Identity(shape, shape);
    this->_shape = shape;
}

GLR GLR::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{GLR} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The GLR element.
    */

    return GLR(shape);
}

Eigen::MatrixXd GLR::project(const Eigen::MatrixXd& other)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in GLR \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(other.rows(), other.cols());
    return other(Eigen::seqN(0, shape), Eigen::seqN(0, shape));
}

size_t GLR::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    return static_cast<size_t>(std::pow(this->_shape, 2));
}

size_t GLR::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return this->_shape;
}

size_t GLR::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    return static_cast<size_t>(std::pow(this->_shape, 2));
}

Eigen::MatrixXd GLR::get_matrix() const
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

GLR GLR::inverse() const
{
    /*! \f{equation*}{ (GLR) \rightarrow GLR \f}
    * 
    * Returns the inverse.
    */

    return this->data.inverse();
}

// Data representation

Eigen::VectorXd GLR::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */

    return this->data.reshaped<Eigen::RowMajor>();
}

void GLR::unserialize(const Eigen::VectorXd& vec)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the GLR object from a serialized vector.
    */

    const size_t vdim = vec.size();
    const size_t max_ind = std::min(this->get_dimension(), vdim);

    for (size_t vind = 0; vind < max_ind; vind++)
    {
        const size_t row = static_cast<size_t>(std::floor(vind / this->_shape));
        const size_t col = vind % this->_shape;
        this->data(row, col) = vec(vind);
    }
}

void GLR::unserialize(std::initializer_list<double> vec)
{
    this->unserialize(Eigen::VectorXd{std::move(vec)});
}

double GLR::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
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

GLR GLR::operator*(const GLR & other) const
{
    /*! \f{equation*}{ (GLR, GLR) \rightarrow GLR \f}
    *
    * Group product.
    */

    assert(this->_shape == other.get_shape());
    return GLR(this->data * other.data);
}

GLR & GLR::operator*=(const GLR & other)
{
    /*! \f{equation*}{ (GLR, GLR) \rightarrow GLR \f}
    *
    * In place group product.
    */

    assert(this->_shape == other.get_shape());
    this->data *= other.data;
    return *this;
}

std::ostream & operator<<(std::ostream& os, const GLR & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << static_cast<const Eigen::MatrixXd>(other.data);
    return os;
}

}
