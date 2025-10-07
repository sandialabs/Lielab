#include "rn.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string rn::to_string() const
{
    if (this->_shape == 0) return "r^nan";
    return "r^" + std::to_string(this->_shape-1);
}

rn::rn() : rn(0)
{
    /*! \f{equation*}{() \rightarrow \mathfrak{rn} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::rn x, y, z;
    * 
    */

}

rn::rn(const size_t n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{rn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{rn}\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::rn x(3), y(4), z(5);
    * 
    * @param[in] shape The shape of the data matrix.
    */

    this->_shape = n + 1;
    this->data = Eigen::VectorXd::Zero(n);
}

rn rn::basis(const ptrdiff_t i, const size_t n)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{rn} \f}
    *
    * Returns the i'th basis element of the rn algebra.
    * 
    * @param[in] i The basis vector.
    * @param[in] n The size of the algebra.
    * @param[out] out The rn element. 
    */

    rn out(n);
    if (i < 0) return out;

    const size_t ind = static_cast<size_t>(i);

    const size_t dim = out.get_dimension();
    if (ind >= dim) return out;

    out.data(ind) = 1.0;

    return out;
}

rn rn::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{rn} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The rn element. 
    */

    if (shape == 0)
    {
        rn out;
        out._shape = 0;
        return out;
    }

    return rn(shape - 1);
}

size_t rn::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    if (this->_shape == 0) return 0; // Return nan?

    return this->_shape - 1;
}

Eigen::VectorXd rn::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */

    return this->data;
}

void rn::set_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ \mathfrak{rn} := \mathbb{R}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    */

    const size_t vdim = vector.size();
    const size_t max_ind = std::min(this->get_dimension(), vdim);

    for (size_t vind = 0; vind < max_ind; vind++)
    {
        this->data(vind) = vector(vind);
    }
}

void rn::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

rn::matrix_t rn::get_matrix() const
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

    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(this->_shape, this->_shape);

    if (this->_shape == 0) return out;

    for (size_t ii = 0; ii < this->_shape-1; ii++)
    {
        out(ii, this->_shape - 1) = data(ii);
    }

    return out;
}

double rn::operator()(const ptrdiff_t index) const
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

    return this->data(_index);
}

double rn::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
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

    if (_index1 == shape - 1) return 0.0;
    if (_index2 != shape - 1) return 0.0;

    return this->data(_index1);
}

rn rn::operator+(const rn & other) const
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape - 1);
    const Eigen::VectorXd new_vector = this->data(slice) + other.data(slice);
    return rn::from_vector(new_vector);
}

rn & rn::operator+=(const rn & other)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    assert(this->_shape == other.get_shape());
    this->data += other.data;
    return *this;
}

rn rn::operator-(const rn & other) const
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape - 1);
    const Eigen::VectorXd new_vector = this->data(slice) - other.data(slice);
    return rn::from_vector(new_vector);
}

rn & rn::operator-=(const rn & other)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    assert(this->_shape == other.get_shape());
    this->data -= other.data;
    return *this;
}

rn rn::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * Unary negative of the vector.
    */

    return rn::from_vector(-this->data);
}

rn rn::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
    *
    * Scalar product.
    */

    return rn::from_vector(this->data * other);
}

rn operator*(const double other, const rn & rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * Scalar product.
    */

    return rn::from_vector(rhs.data * other);
}

rn & rn::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
    *
    * In place scalar product.
    */

    this->data *= other;
    return *this;
}

rn rn::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
    *
    * Scalar division.
    */

    return rn::from_vector(this->data / other);
}

rn & rn::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
    *
    * In place scalar division.
    */

    this->data /= other;
    return *this;
}

rn rn::from_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{rn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{rn}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    const size_t shape = vector.size() + 1;
    rn out = rn::from_shape(shape);
    out.set_vector(vector);

    return out;
}

rn rn::from_vector(std::initializer_list<double> other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{rn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{rn}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    return rn::from_vector(Eigen::VectorXd{std::move(other)});
}

Eigen::MatrixXd rn::project(const Eigen::MatrixXd & other)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{rn} \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(other.rows(), other.cols());

    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(shape, shape);

    for (size_t ii = 0; ii < shape - 1; ii++)
    {
        out(ii, shape-1) = other(ii, shape-1);
    }

    return out;
}

std::ostream& operator<<(std::ostream & os, const rn & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    
    os << static_cast<const Eigen::VectorXd>(other.data);
    return os;
}

}
