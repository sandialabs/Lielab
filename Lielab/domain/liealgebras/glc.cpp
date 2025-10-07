#include "glc.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string glc::to_string() const
{
    return "gl(" + std::to_string(this->_shape) + ", C)";
}

glc::glc() : glc(0)
{
    /*! \f{equation*}{() \rightarrow \mathfrak{glc} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::glc x, y, z;
    * 
    */

}

glc::glc(const size_t n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{glc} \f}
    *
    * Constructor instantiating an \f$\mathfrak{glc}\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::glc x(3), y(4), z(5);
    * 
    * @param[in] n The shape of the data matrix.
    */

    this->data = Eigen::MatrixXcd::Zero(n, n);
    this->_shape = n;
}

glc glc::basis(const ptrdiff_t i, const size_t n) 
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{glc} \f}
    *
    * Returns the i'th basis element of the glc algebra.
    * 
    * @param[in] i The basis vector.
    * @param[in] n The shape of the algebra.
    * @param[out] out The glc element.
    */

    glc out(n);
    if (i < 0) return out;

    const size_t ind = static_cast<size_t>(i);

    if (ind >= 2*n*n) return out;

    const size_t row = static_cast<size_t>(std::floor((ind/2) / n));
    const size_t col = (ind/2) % n;

    if ((ind % 2) == 0)
    {
        out.data(row, col).real(1.0);
    }
    else
    {
        out.data(row, col).imag(1.0);
    }

    return out;
}

glc glc::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{glc} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The cn element. 
    */

    return glc(shape);
}

size_t glc::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    return static_cast<size_t>(2*std::pow(this->_shape, 2));
}

size_t glc::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the algebra.
    */

    return this->_shape;
}

Eigen::VectorXd glc::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */
    
    const size_t dim = this->get_dimension();
    const Eigen::MatrixXcd A = this->get_matrix();

    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);
    size_t kk = 0;
    for (size_t ii = 0; ii < this->_shape; ii++)
    {
        for (size_t jj = 0; jj < this->_shape; jj++)
        {
            out(kk) = std::real(A(ii,jj));
            out(kk+1) = std::imag(A(ii,jj));
            kk = kk + 2;
        }
    }

    return out;
}

void glc::set_vector(const Eigen::VectorXd& vector) 
{
    /*! \f{equation*}{ \mathfrak{glc} := \mathbb{R}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    */

    const size_t vdim = vector.size();
    const size_t max_ind = std::min(this->get_dimension(), vdim);
    
    for (size_t vind = 0; vind < max_ind; vind++)
    {
        const size_t rem = vind % 2;
        const size_t row = static_cast<size_t>(std::floor((vind/2) / this->_shape));
        const size_t col = (vind/2) % this->_shape;
        
        if (rem == 0)
        {
            this->data(row, col).real(vector(vind));
        }
        else
        {
            this->data(row, col).imag(vector(vind));
        }
    }
}

void glc::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
    
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

glc::matrix_t glc::get_matrix() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times n} \f}
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

    return this->data.cast<std::complex<double>>();
}

double glc::operator()(const ptrdiff_t index) const
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
    
    const size_t row = static_cast<size_t>(std::floor((_index/2) / shape));
    const size_t col = (_index/2) % shape;

    if ((_index % 2) == 0)
    {
        return this->data(row, col).real();
    }

    return this->data(row, col).imag();
}

std::complex<double> glc::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    * Gets a value in the square matrix representation.
    */
    
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t shape = this->get_shape();
    if (shape == 0) return std::complex<double>(nan, nan);

    if (index1 >= static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);
    if (index2 >= static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);
    if (std::abs(index1) > static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);
    if (std::abs(index2) > static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);

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

glc glc::operator+(const glc & other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = this->get_matrix();
    const Eigen::MatrixXcd rhs_matrix = other.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return glc(new_matrix);
}

glc & glc::operator+=(const glc & other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd other_matrix = other.get_matrix();

    this->data += other_matrix(slice, slice);
    return *this;
}

glc glc::operator-(const glc & other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = this->get_matrix();
    const Eigen::MatrixXcd rhs_matrix = other.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return glc(new_matrix);
}

glc & glc::operator-=(const glc & other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd other_matrix = other.get_matrix();

    this->data -= other_matrix(slice, slice);
    return *this;
}

glc glc::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Unary negative of the vector.
    */

    const Eigen::MatrixXcd this_matrix = this->get_matrix();

    return -this_matrix;
}

glc glc::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{R}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar product.
    */

    const Eigen::MatrixXcd this_matrix = this->get_matrix();

    return this_matrix * other;
}

glc operator*(const double other, const glc & rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

glc glc::operator*(const std::complex<double> other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar product.
    */

    const Eigen::MatrixXcd this_matrix = this->get_matrix();

    return this_matrix * other;
}

glc operator*(const std::complex<double> other, const glc & rhs)
{
    /*! \f{equation*}{ (\mathbb{C}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

glc & glc::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{R}) \rightarrow \mathfrak{glc} \f}
    *
    * In place scalar product.
    */

    this->data *= other;
    return *this;
}

glc & glc::operator*=(const std::complex<double> other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
    *
    * In place scalar product.
    */

    this->data *= other;
    return *this;
}

glc glc::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{R}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar division.
    */

    const Eigen::MatrixXcd this_matrix = this->get_matrix();
    return this_matrix / other;
}

glc glc::operator/(const std::complex<double> other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar division.
    */

    const Eigen::MatrixXcd this_matrix = this->get_matrix();
    return this_matrix / other;
}

glc & glc::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{R}) \rightarrow \mathfrak{glc} \f}
    *
    * In place scalar division.
    */

    this->data /= other;
    return *this;
}

glc & glc::operator/=(const std::complex<double> other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
    *
    * In place scalar division.
    */

    this->data /= other;
    return *this;
}

glc glc::from_vector(const Eigen::VectorXd& other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{glc} \f}
    *
    * Constructor instantiating a \f$\mathfrak{glc}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    const size_t len = other.size();
    const size_t shape = static_cast<size_t>(std::ceil(std::sqrt(static_cast<double>(len)/2.0)));
    glc out(shape);
    out.set_vector(other);
    return out;
}

glc glc::from_vector(std::initializer_list<double> other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{glc} \f}
    *
    * Constructor instantiating an \f$\mathfrak{glc}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    return glc::from_vector(Eigen::VectorXd{std::move(other)});
}

glc glc::from_complex_vector(const Eigen::VectorXcd &other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{glc} \f}
    *
    * Constructor instantiating a \f$\mathfrak{glc}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    const size_t len = other.size();

    const double shape = std::ceil(std::sqrt(static_cast<double>(len)));

    Eigen::VectorXd temp = Eigen::VectorXd::Zero(static_cast<size_t>(2*shape*shape));

    for (size_t ii = 0; ii < len; ii++)
    {
        temp(2*ii) = other(ii).real();
        temp(2*ii+1) = other(ii).imag();
    }

    glc out(static_cast<size_t>(shape));
    out.set_vector(temp);

    return out;
}

glc glc::from_complex_vector(std::initializer_list<std::complex<double>> other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{glc} \f}
    *
    * Constructor instantiating an \f$\mathfrak{glc}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    return glc::from_complex_vector(Eigen::VectorXcd{std::move(other)});
}

Eigen::MatrixXcd glc::project(const Eigen::MatrixXcd & other)
{
    /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in \mathfrak{glc} \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(other.rows(), other.cols());
    return other(Eigen::seqN(0, shape), Eigen::seqN(0, shape));
}

std::ostream& operator<<(std::ostream & os, const glc & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    
    os << static_cast<const Eigen::MatrixXcd>(other.data);
    return os;
}

}
