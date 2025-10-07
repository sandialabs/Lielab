#include "glr.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string glr::to_string() const
{
    return "gl(" + std::to_string(this->_shape) + ", R)";
}

glr::glr() : glr(0)
{
    /*! \f{equation*}{() \rightarrow \mathfrak{glr} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::glr x, y, z;
    * 
    */

}

glr::glr(const size_t n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{glr} \f}
    *
    * Constructor instantiating an \f$\mathfrak{glr}\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::glr x(3), y(4), z(5);
    * 
    * @param[in] n The shape of the data matrix.
    */

    this->data = Eigen::MatrixXd::Zero(n, n);
    this->_shape = n;
}

glr glr::basis(const ptrdiff_t i, const size_t n) 
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{glr} \f}
    *
    * Returns the i'th basis element of the glr algebra.
    * 
    * @param[in] i The basis vector.
    * @param[in] n The shape of the algebra.
    * @param[out] out The glr element.
    */

    glr out(n);
    if (i < 0) return out;

    const size_t ind = static_cast<size_t>(i);

    if (ind >= n*n) return out;

    const size_t row = static_cast<size_t>(std::floor(ind / n));
    const size_t col = ind % n;

    out.data(row, col) = 1.0;

    return out;
}

glr glr::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{glr} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The cn element. 
    */

    return glr(shape);
}

size_t glr::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    return static_cast<size_t>(std::pow(this->_shape, 2));
}

size_t glr::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the algebra.
    */

    return this->_shape;
}

Eigen::VectorXd glr::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */
    
    return this->data.reshaped<Eigen::RowMajor>();
}

void glr::set_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ \mathfrak{glr} := \mathbb{R}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    */

    const size_t vdim = vector.size();
    const size_t max_ind = std::min(this->get_dimension(), vdim);

    for (size_t vind = 0; vind < max_ind; vind++)
    {
        const size_t row = static_cast<size_t>(std::floor(vind / this->_shape));
        const size_t col = vind % this->_shape;
        this->data(row, col) = vector(vind);
    }
}

void glr::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

glr::matrix_t glr::get_matrix() const
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

    return data;
}

double glr::operator()(const ptrdiff_t index) const
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the vector representation.
    */
    
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t dim = this->get_dimension();

    if (index >= static_cast<ptrdiff_t>(dim)) return nan;
    if (std::abs(index) > static_cast<ptrdiff_t>(dim)) return nan;

    size_t _index;
    if (index < 0)
    {
        _index = static_cast<size_t>(static_cast<ptrdiff_t>(dim) + index);
    }
    else
    {
        _index = static_cast<size_t>(index);
    }

    const size_t shape = this->get_shape();

    const size_t row = static_cast<size_t>(std::floor(_index / shape));
    const size_t col = _index % shape;

    return this->data(row, col);
}

double glr::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
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

glr glr::operator+(const glr & other) const
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other._shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = this->get_matrix();
    const Eigen::MatrixXd rhs_matrix = other.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return glr(new_matrix);
}

glr & glr::operator+=(const glr & other)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    assert(this->_shape == other.get_shape());
    this->data += other.data;
    return *this;
}

glr glr::operator-(const glr & other) const
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = this->get_matrix();
    const Eigen::MatrixXd rhs_matrix = other.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return glr(new_matrix);
}

glr & glr::operator-=(const glr & other)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    assert(this->_shape == other.get_shape());
    this->data -= other.data;
    return *this;
}

glr glr::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * Unary negative of the vector.
    */

    return -this->data;
}

glr glr::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    *
    * Scalar product.
    */

    return this->data * other;
}

glr operator*(const double other, const glr & rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

glr & glr::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    *
    * In place scalar product.
    */

    this->data *= other;
    return *this;
}

glr glr::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    *
    * Scalar division.
    */

    return this->data / other;
}

glr & glr::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    *
    * In place scalar division.
    */

    this->data /= other;
    return *this;
}

glr glr::from_vector(const Eigen::VectorXd &other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{glr} \f}
    *
    * Constructor instantiating a \f$\mathfrak{glr}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    const size_t len = other.size();
    const size_t shape = static_cast<size_t>(std::ceil(std::sqrt(static_cast<double>(len))));
    glr out(shape);
    out.set_vector(other);
    return out;
}

glr glr::from_vector(std::initializer_list<double> other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{glr} \f}
    *
    * Constructor instantiating an \f$\mathfrak{glr}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    return glr::from_vector(Eigen::VectorXd{std::move(other)});
}

Eigen::MatrixXd glr::project(const Eigen::MatrixXd& other)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{glr} \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(other.rows(), other.cols());
    return other(Eigen::seqN(0, shape), Eigen::seqN(0, shape));
}


std::ostream& operator<<(std::ostream & os, const glr & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    
    os << static_cast<const Eigen::MatrixXd>(other.data);
    return os;
}

}
