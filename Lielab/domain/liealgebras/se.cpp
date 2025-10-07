#include "se.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string se::to_string() const
{
    if (this->_shape == 0) return "se(nan)";
    return "se(" + std::to_string(this->_shape-1) + ")";
}

se::se() : se(0)
{
    /*! \f{equation*}{() \rightarrow \mathfrak{se} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::se x, y, z;
    * 
    */

}

se::se(const size_t n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{se} \f}
    *
    * Constructor instantiating an \f$\mathfrak{se}\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::se x(3), y(4), z(5);
    * 
    * @param[in] shape The shape of the data matrix.
    */
    
    this->_shape = n + 1;
    this->data = Eigen::MatrixXd::Zero(n + 1, n + 1);
}

se se::basis(const ptrdiff_t i, const size_t n)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{se} \f}
    *
    * Returns the i'th basis element of the se algebra.
    * 
    * @param[in] i The basis vector.
    * @param[in] n The size of the algebra.
    * @param[out] out The se element.
    */

    se out(n);
    if (i < 0) return out;

    const size_t ind = static_cast<size_t>(i);

    const size_t dim = out.get_dimension();
    if (ind >= dim) return out;

    Eigen::VectorXd v = Eigen::VectorXd::Zero(dim);
    v(ind) = 1.0;
    out.set_vector(v);
    return out;
}

se se::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{se} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The se element. 
    */

    if (shape == 0)
    {
        se out;
        out._shape = 0;
        out.data = Eigen::MatrixXd::Zero(0,0);
        return out;
    }

    return se(shape - 1);
}

size_t se::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    if (this->_shape == 0) return 0; // Return nan?
    if (this->_shape == 1) return 0; // Return nan?

    return (this->_shape - 1) * (this->_shape - 2) / 2 + this->_shape - 1;
}

Eigen::VectorXd se::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */

    const size_t dim = this->get_dimension();
    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);
    if (this->_shape == 0) return out;

    int k = 0;

    for (size_t ii = 0; ii < this->_shape - 1; ii++)
    {
        out(k) = this->data(ii, this->_shape - 1);
        k++;
    }

    if (dim <= 1)
    {
        return out;
    }

    for (size_t ii = this->_shape - 2; ii > 0; ii--)
    {
        for (size_t jj = this->_shape - 1; jj > ii; jj--)
        {
            out(k) = this->data(ii-1, jj-1)/std::pow(-1.0, ii+jj);
            k++;
        }
    }

    return out;
}

void se::set_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ \mathfrak{se} := \mathbb{R}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    */

    const size_t vdim = vector.size();
    const size_t max_ind = std::min(this->get_dimension(), vdim);

    size_t k = 0;
    if (k >= max_ind) return;

    for (size_t ii = 0; ii < this->_shape - 1; ii++)
    {
        this->data(ii, this->_shape - 1) = vector(k);
        k++;
        if (k >= max_ind) return;
    }

    for (size_t ii = this->_shape - 2; ii > 0; ii--)
    {
        for (size_t jj = this->_shape - 1; jj > ii; jj--)
        {
            this->data(ii-1, jj-1) =  std::pow(-1, (ii+jj))*vector(k);
            this->data(jj-1, ii-1) = -std::pow(-1, (ii+jj))*vector(k);
            k++;
            if (k >= max_ind) return;
        }
    }
}

void se::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

se::matrix_t se::get_matrix() const
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

double se::operator()(const ptrdiff_t index) const
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the column vector representation.
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

    const Eigen::VectorXd vector = this->get_vector();
    return vector(_index);
}

double se::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t shape = this->get_shape();
    if (shape <= 1) return nan;

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

se se::operator+(const se & other) const
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = this->get_matrix();
    const Eigen::MatrixXd rhs_matrix = other.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return se(new_matrix);
}

se & se::operator+=(const se & other)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    assert(this->_shape == other.get_shape());
    this->data += other.data;
    return *this;
}

se se::operator-(const se & other) const
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = this->get_matrix();
    const Eigen::MatrixXd rhs_matrix = other.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return se(new_matrix);
}

se & se::operator-=(const se & other)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    assert(this->_shape == other.get_shape());
    this->data -= other.data;
    return *this;
}

se se::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * Unary negative of the vector.
    */
    
    return -this->data;
}

se se::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
    *
    * Scalar product.
    */

    se out = this->data * other;
    return out;
}

se operator*(const double other, const se & rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

se & se::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
    *
    * In place scalar multiplication.
    */

    this->data *= other;
    return *this;
}

se se::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
    *
    * Scalar division.
    */

    Eigen::MatrixXd out = this->data / other;
    return out;
}

se & se::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{RR}) \rightarrow \mathfrak{se} \f}
    *
    * In place scalar division.
    */

    this->data /= other;
    return *this;
}

se se::from_vector(const Eigen::VectorXd& other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{se} \f}
    *
    * Constructor instantiating an \f$\mathfrak{se}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    const size_t len = other.size();
    const size_t shape = static_cast<size_t>(std::ceil(std::sqrt(8.0*len + 1.0)/2.0 + 0.5));

    se out = se::from_shape(shape);
    out.set_vector(other);
    return out;
}

se se::from_vector(std::initializer_list<double> other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{se} \f}
    *
    * Constructor instantiating an \f$\mathfrak{se}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    return se::from_vector(Eigen::VectorXd{std::move(other)});
}

// static Eigen::MatrixXd project(const Eigen::MatrixXd & other)
// {
//     /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{se} \f}
//     *
//     * Projects a matrix suitable for data.
//     */

//     if (other.rows() != other.cols())
//     {
//         throw Lielab::utils::Error("Size of the matrix must be square.");
//     }

//     return (other - other.transpose())/2;
// }

// Eigen::MatrixXd se::project(const Eigen::MatrixXd & other)
// {
//     /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{se} \f}
//     *
//     * Projects a matrix suitable for data.
//     */

//     if (other.rows() != other.cols())
//     {
//         throw Lielab::utils::Error("Size of the matrix must be square.");
//     }

//     return (other - other.transpose())/2;
// }

std::ostream & operator<<(std::ostream & os, const se & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    
    os << static_cast<const Eigen::MatrixXd>(other.data);
    return os;
}

}

