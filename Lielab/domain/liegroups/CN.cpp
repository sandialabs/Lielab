#include "CN.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string CN::to_string() const
{
    const size_t shape = this->get_shape();
    if (shape == 0) return "C^nan";
    return "C^" + std::to_string(shape-1);
}

CN::CN() : CN(0)
{
    /*! \f{equation*}{ () \rightarrow CN \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::CN x, y, z;
    * 
    */

}

CN::CN(const size_t n)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow CN \f}
    *
    * Constructor instantiating an \f$CN\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::CN x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */
    
    this->_shape = n + 1;
    this->data = Eigen::VectorXcd::Zero(n);
}

CN CN::from_shape(const size_t shape)
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
        CN out;
        out._shape = 0;
        return out;
    }

    return CN(shape-1);
}

size_t CN::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    if (this->_shape == 0) return 0; // TODO: Return nan?
    return 2*(this->_shape - 1);
}

size_t CN::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return this->_shape;
}

size_t CN::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    if (this->_shape == 0) return 0;
    return 2*(this->_shape - 1);
}

Eigen::MatrixXcd CN::get_matrix() const
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

    Eigen::MatrixXcd out = Eigen::MatrixXcd::Identity(this->_shape, this->_shape);

    if (this->_shape == 0) return out;

    for (size_t ii = 0; ii < this->_shape - 1; ii++)
    {
        out(ii, this->_shape - 1) = this->data(ii);
    }

    return out;
}

CN CN::inverse() const
{
    /*! \f{equation*}{ (CN) \rightarrow CN \f}
    * 
    * Returns the inverse.
    */

    return CN::from_complex_vector(-this->data);
}

// Data representation

Eigen::VectorXd CN::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
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

void CN::unserialize(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ (\mathbb{C}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the CN object from a serialized vector.
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

void CN::unserialize(const std::initializer_list<double> other)
{
    /*!
    */

    this->unserialize(Eigen::VectorXd{std::move(other)});
}

double CN::operator()(const ptrdiff_t index) const
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

// std::complex<double> & CN::operator()(const size_t index)
// {
//     /*! \f{equation*}{ CN(\mathbb{Z}) := \mathbb{C} \f}
//     *
//     * Assignment of a value in the column vector representation.
//     */

//     return this->data(index);
// }

std::complex<double> CN::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
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

    if (_index1 == _index2) return std::complex<double>(1.0, 0.0);

    if (_index1 == shape - 1) return std::complex<double>(0.0, 0.0);
    if (_index2 != shape - 1) return std::complex<double>(0.0, 0.0);

    return this->data(_index1);
}

std::complex<double> CN::operator[](const ptrdiff_t index) const
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

CN CN::operator*(const CN & other) const
{
    /*! \f{equation*}{ (CN, CN) \rightarrow CN \f}
    *
    * Group product.
    */

    assert(this->_shape == other.get_shape());

    return CN::from_complex_vector(this->data + other.data);
}

CN & CN::operator*=(const CN & other)
{
    /*! \f{equation*}{ (CN, CN) \rightarrow CN \f}
    *
    * In place group product.
    */

    assert(this->_shape == other.get_shape());

    this->data += other.data;
    return *this;
}

CN CN::from_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{CN} \f}
    *
    * Constructor instantiating an \f$\mathfrak{CN}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] vector The object to instantiate from as an imaginary vector.
    */

    const size_t shape = static_cast<size_t>(std::ceil(vector.size()/2.0)) + 1;
    CN out = CN::from_shape(shape);
    out.unserialize(vector);

    return out;
}

CN CN::from_vector(std::initializer_list<double> vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{CN} \f}
    *
    * Constructor instantiating an \f$\mathfrak{CN}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    return CN::from_vector(Eigen::VectorXd{std::move(vector)});
}

CN CN::from_complex_vector(const Eigen::VectorXcd &other)
{
    /*! \f{equation}{(\mathbb{C}^{n \times 1}) \rightarrow \mathfrak{CN} \f}
    *
    * Constructor instantiating an \f$\mathfrak{CN}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    const size_t n = other.size();
    CN out(n);
    out.data = other;
    return out;
}

CN CN::from_complex_vector(const std::initializer_list<std::complex<double>> other)
{
    /*! \f{equation}{(\mathbb{C}^{n \times 1}) \rightarrow \mathfrak{CN} \f}
    *
    * Constructor instantiating an \f$\mathfrak{CN}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    return CN::from_complex_vector(Eigen::VectorXcd{std::move(other)});
}

Eigen::VectorXcd CN::to_complex_vector() const
{
    /*! \f{equation}{ () \rightarrow (\mathbb{C}^{n \times 1}) \f}
    *
    * @param[out] vec The object as a complex vector.
    */

    return Eigen::VectorXcd(this->data);
}

Eigen::MatrixXcd CN::project(const Eigen::MatrixXcd & other)
{
    /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in CN \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(other.rows(), other.cols());

    Eigen::MatrixXcd out = Eigen::MatrixXcd::Identity(shape, shape);

    for (size_t ii = 0; ii < shape - 1; ii++)
    {
        out(ii, shape-1) = other(ii, shape-1);
    }

    return out;
}

std::ostream & operator<<(std::ostream& os, const CN & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << static_cast<const Eigen::VectorXcd>(other.data);
    return os;
}

}
