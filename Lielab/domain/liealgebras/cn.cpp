#include "cn.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <complex>

namespace Lielab::domain
{

bool cn::is_abelian() const
{
    return true;
}

std::string cn::to_string() const
{
    const size_t shape = this->get_shape();
    if (shape == 0) return "c^nan";
    return "c^" + std::to_string(shape-1);
}

cn::cn() : cn(0)
{
    /*! \f{equation*}{() \rightarrow \mathfrak{cn} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::cn x, y, z;
    * 
    */

}

cn::cn(const size_t n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{cn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{cn}\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::cn x(3), y(4), z(5);
    * 
    * @param[in] shape The shape of the data matrix.
    */

    this->_shape = n + 1;
    this->data = Eigen::VectorXcd::Zero(n);
}

cn cn::basis(const ptrdiff_t i, const size_t n)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{cn} \f}
    *
    * Returns the i'th basis element of the cn algebra.
    * 
    * @param[in] i The basis vector.
    * @param[in] n The size of the algebra.
    * @param[out] out The cn element.
    */

    cn out(n);
    if (i < 0) return out;

    const size_t dim = out.get_dimension();
    const size_t ind = static_cast<size_t>(i);

    if (ind >= dim) return out;

    const size_t indz = ind/2;
    const size_t rem = ind%2;
    if (rem == 0) out.data(indz) = std::complex<double>(1.0, 0.0);
    if (rem == 1) out.data(indz) = std::complex<double>(0.0, 1.0);

    return out;
}

cn cn::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{cn} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The cn element. 
    */

    if (shape == 0)
    {
        cn out;
        out._shape = 0;
        return out;
    }

    return cn(shape-1);
}

size_t cn::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    if (this->_shape == 0) return 0; // Return nan??

    return 2*(this->_shape - 1);
}

Eigen::VectorXd cn::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */

    const size_t dim = this->get_dimension();
    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);

    size_t kk = 0;
    for (size_t ii = 0; ii < dim/2; ii++)
    {
        out(kk) = std::real(this->data(ii));
        kk += 1;
        out(kk) = std::imag(this->data(ii));
        kk += 1;
    }

    return out;
}

void cn::set_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ \mathfrak{cn} := \mathbb{C}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    */

    const size_t vdim = vector.size();
    const size_t max_ind = std::min(this->get_dimension(), vdim);

    size_t dind = 0;
    for (size_t vind = 0; vind < max_ind; vind++)
    {
        const size_t rem = vind % 2;

        if (rem == 0)
        {
            this->data(dind).real(vector(vind));
        }
        else
        {
            this->data(dind).imag(vector(vind));
            dind += 1;
        }
    }
}

void cn::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

cn::matrix_t cn::get_matrix() const
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

    Eigen::MatrixXcd out = Eigen::MatrixXcd::Zero(this->_shape, this->_shape);

    if (this->_shape == 0) return out;

    for (size_t ii = 0; ii < this->_shape - 1; ii++)
    {
        out(ii, this->_shape - 1) = this->data(ii);
    }

    return out;
}

double cn::operator()(const ptrdiff_t index) const
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

    if (_index % 2 == 0) return this->data(_index/2).real();
    return this->data(_index/2).imag();
}

std::complex<double> cn::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
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

    if (_index1 == shape - 1) return std::complex<double>(0.0, 0.0);
    if (_index2 != shape - 1) return std::complex<double>(0.0, 0.0);

    return this->data(_index1);
}

std::complex<double> cn::operator[](const ptrdiff_t index) const
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t dim = this->get_dimension();

    if (index >= static_cast<ptrdiff_t>(dim/2)) return std::complex<double>(nan, nan);
    if (std::abs(index) > static_cast<ptrdiff_t>(dim/2)) return std::complex<double>(nan, nan);

    size_t _index;
    if (index < 0)
    {
        _index = static_cast<size_t>(static_cast<ptrdiff_t>(dim/2) + index);
    }
    else
    {
        _index = static_cast<size_t>(index);
    }

    return this->data(_index);
}

cn cn::operator+(const cn & other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape - 1);
    const Eigen::VectorXcd new_vector = this->data(slice) + other.data(slice);
    return cn::from_complex_vector(new_vector);
}

cn & cn::operator+=(const cn & other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    this->data += other.data;
    return *this;
}

cn cn::operator-(const cn & other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape - 1);
    const Eigen::VectorXcd new_vector = this->data(slice) - other.data(slice);
    return cn::from_complex_vector(new_vector);
}

cn & cn::operator-=(const cn & other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    this->data -= other.data;
    return *this;
}

cn cn::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Unary negative of the vector.
    */

    return cn::from_complex_vector(-this->data);
}

cn cn::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */

    return cn::from_complex_vector(this->data*other);
}

cn operator*(const double other, const cn & rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */

    return cn::from_complex_vector(other*rhs.data);
}

cn cn::operator*(const std::complex<int> other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */
    
    const std::complex<double> otherd = std::complex<double>(static_cast<double>(other.real()),
                                                             static_cast<double>(other.imag()));
    return cn::from_complex_vector(this->data*otherd);
}

cn cn::operator*(const std::complex<double> other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */

    return cn::from_complex_vector(this->data*other);
}

cn operator*(const std::complex<int> other, const cn & rhs)
{
    /*! \f{equation*}{ (\mathbb{C}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */
    
    const std::complex<double> otherd = std::complex<double>(static_cast<double>(other.real()),
                                                             static_cast<double>(other.imag()));
    return cn::from_complex_vector(otherd*rhs.data);
}

cn operator*(const std::complex<double> other, const cn & rhs)
{
    /*! \f{equation*}{ (\mathbb{C}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */

    return cn::from_complex_vector(other*rhs.data);
}

cn & cn::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar product.
    */

    this->data *= other;
    return *this;
}

cn & cn::operator*=(const std::complex<int> other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar product.
    */
    
    const std::complex<double> otherd = std::complex<double>(static_cast<double>(other.real()),
                                                             static_cast<double>(other.imag()));
    this->data *= otherd;
    return *this;
}

cn & cn::operator*=(const std::complex<double> other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar product.
    */

    this->data *= other;
    return *this;
}

cn cn::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar division.
    */

    return cn::from_complex_vector(this->data/other);
}

cn cn::operator/(const std::complex<int> other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar division.
    */

    const std::complex<double> otherd = std::complex<double>(static_cast<double>(other.real()),
                                                             static_cast<double>(other.imag()));
    return cn::from_complex_vector(this->data/otherd);
}

cn cn::operator/(const std::complex<double> other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar division.
    */

    return cn::from_complex_vector(this->data/other);
}

cn & cn::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar division.
    */

    this->data /= other;
    return *this;
}

cn & cn::operator/=(const std::complex<int> other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar division.
    */

    const std::complex<double> otherd = std::complex<double>(static_cast<double>(other.real()),
                                                             static_cast<double>(other.imag()));
    this->data /= otherd;
    return *this;
}

cn & cn::operator/=(const std::complex<double> other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar division.
    */

    this->data /= other;
    return *this;
}

cn cn::from_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{cn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{cn}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] vector The object to instantiate from as an imaginary vector.
    */

    const size_t shape = static_cast<size_t>(std::ceil(vector.size()/2.0)) + 1;
    cn out = cn::from_shape(shape);
    out.set_vector(vector);

    return out;
}

cn cn::from_vector(std::initializer_list<double> vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{cn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{cn}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    return cn::from_vector(Eigen::VectorXd{std::move(vector)});
}

cn cn::from_complex_vector(const Eigen::VectorXcd &other)
{
    /*! \f{equation}{(\mathbb{C}^{n \times 1}) \rightarrow \mathfrak{cn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{cn}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    const size_t n = other.size();
    cn out(n);
    out.data = other;
    return out;
}

cn cn::from_complex_vector(std::initializer_list<std::complex<double>> other)
{
    /*! \f{equation}{(\mathbb{C}^{n \times 1}) \rightarrow \mathfrak{cn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{cn}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    return cn::from_complex_vector(Eigen::VectorXcd{std::move(other)});
}

Eigen::VectorXcd cn::to_complex_vector() const
{
    /*
     * Returns complex vector.
     */
    
    return this->data;
}

Eigen::MatrixXcd cn::project(const Eigen::MatrixXcd & other)
{
    /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in \mathfrak{cn} \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(other.rows(), other.cols());

    Eigen::MatrixXcd out = Eigen::MatrixXcd::Zero(shape, shape);

    for (size_t ii = 0; ii < shape - 1; ii++)
    {
        out(ii, shape-1) = other(ii, shape-1);
    }

    return out;
}

std::ostream& operator<<(std::ostream & os, const cn & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    
    os << static_cast<const Eigen::VectorXcd>(other.data);
    return os;
}

}
