#include "sp.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string sp::to_string() const
{
    return "sp(" + std::to_string(this->_shape) + ", R)";
}

sp::sp() : sp(0)
{
    /*! \f{equation*}{ () \rightarrow \mathfrak{sp} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::sp x, y, z;
    * 
    */

}

sp::sp(const size_t n)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{sp} \f}
    *
    * Constructor instantiating an \f$\mathfrak{sp}\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::sp x(2), y(4), z(6);
    * 
    * @param[in] shape The shape of the data matrix.
    */

    this->_shape = 2*static_cast<size_t>(std::floor(n/2));
    this->data = Eigen::MatrixXd::Zero(this->_shape, this->_shape);
}

sp sp::basis(const ptrdiff_t i, const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{sp} \f}
    *
    * Returns the i'th basis element of the sp algebra.
    * 
    * @param[in] i The basis vector.
    * @param[in] shape The shape of the algebra.
    * @param[out] out The rn element.
    */

    sp out(shape);
    if (i < 0) return out;

    const size_t ind = static_cast<size_t>(i);

    const size_t dim = out.get_dimension();
    if (ind >= dim) return out;

    Eigen::VectorXd v = Eigen::VectorXd::Zero(dim);
    v(ind) = 1.0;
    out.set_vector(v);
    return out;
}

sp sp::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{se} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The se element. 
    */

    return sp(shape);
}

size_t sp::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    return this->_shape * (this->_shape + 1) / 2;
}

Eigen::VectorXd sp::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */

    const size_t dim = this->get_dimension();
    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);

    if (dim == 0) return out;

    size_t k = 0;

    for (size_t ii = 0; ii < this->_shape/2; ii++)
    {
        for (size_t jj = 0; jj < this->_shape/2; jj++)
        {
            out(k) = (this->data(ii, jj) - this->data(this->_shape/2 + jj, this->_shape/2 + ii))/2.0;
            k++;
        }
    }

    for (size_t ii = 0; ii < this->_shape/2; ii++)
    {
        for (size_t jj = ii; jj < this->_shape/2; jj++)
        {
            out(k) = (this->data(ii, this->_shape/2 + jj) + this->data(jj, this->_shape/2 + ii))/2.0;
            k++;
        }
    }

    for (size_t ii = 0; ii < this->_shape/2; ii++)
    {
        for (size_t jj = ii; jj < this->_shape/2; jj++)
        {
            out(k) = (this->data(this->_shape/2 + ii, jj) + this->data(this->_shape/2 + jj, ii))/2.0;
            k++;
        }
    }

    return out;
}

void sp::set_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ \mathfrak{sp} := \mathbb{R}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    * TODO: Set shape.
    */

    const size_t vdim = vector.size();
    const size_t max_ind = std::min(this->get_dimension(), vdim);

    size_t k = 0;
    if (k >= max_ind) return;

    for (size_t ii = 0; ii < this->_shape/2; ii++)
    {
        for (size_t jj = 0; jj < this->_shape/2; jj++)
        {
            this->data(ii, jj) = vector(k);
            this->data(this->_shape/2 + jj, this->_shape/2 + ii) = -vector(k);
            k++;
            if (k >= max_ind) return;
        }
    }

    for (size_t ii = 0; ii < this->_shape/2; ii++)
    {
        for (size_t jj = ii; jj < this->_shape/2; jj++)
        {
            this->data(ii, this->_shape/2 + jj) = vector(k);
            this->data(jj, this->_shape/2 + ii) = vector(k);
            k++;
            if (k >= max_ind) return;
        }
    }

    for (size_t ii = 0; ii < this->_shape/2; ii++)
    {
        for (size_t jj = ii; jj < this->_shape/2; jj++)
        {
            this->data(this->_shape/2 + ii, jj) = vector(k);
            this->data(this->_shape/2 + jj, ii) = vector(k);
            k++;
            if (k >= max_ind) return;
        }
    }
}

void sp::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

sp::matrix_t sp::get_matrix() const
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

double sp::operator()(const ptrdiff_t index) const
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

double sp::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t shape = this->get_shape();
    if (shape == 0) return nan;
    if (shape % 2 == 1) return nan;

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

sp sp::operator+(const sp & other) const
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = this->get_matrix();
    const Eigen::MatrixXd rhs_matrix = other.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return sp(new_matrix);
}

sp & sp::operator+=(const sp & other)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    assert(this->_shape == other.get_shape());
    this->data += other.data;
    return *this;
}

sp sp::operator-(const sp & other) const
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = this->get_matrix();
    const Eigen::MatrixXd rhs_matrix = other.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return sp(new_matrix);
}

sp & sp::operator-=(const sp & other)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    assert(this->_shape == other.get_shape());
    this->data -= other.data;
    return *this;
}

sp sp::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * Unary negative of the vector.
    */

    Eigen::MatrixXd out = -this->data;
    return out;
}

sp sp::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
    *
    * Scalar product.
    */

    Eigen::MatrixXd out = this->data * other;
    return out;
}

sp operator*(const double other, const sp & rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

sp & sp::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
    *
    * In place scalar product.
    */

    this->data *= other;
    return *this;
}

sp sp::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
    *
    * Scalar division.
    */

    Eigen::MatrixXd out = this->data / other;
    return out;
}

sp & sp::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
    *
    * In place scalar division.
    */

    this->data /= other;
    return *this;
}

sp sp::from_vector(const Eigen::VectorXd& other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{sp} \f}
    *
    * Constructor instantiating an \f$\mathfrak{sp}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    const size_t len = other.size();
    size_t shape = static_cast<size_t>(std::ceil((-0.5 + std::sqrt(0.25 + 2*len))/(2*0.5)));
    if (shape % 2 != 0) shape++;

    sp out(shape);
    out.set_vector(other);
    return out;
}

sp sp::from_vector(std::initializer_list<double> other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{sp} \f}
    *
    * Constructor instantiating an \f$\mathfrak{sp}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    return sp::from_vector(Eigen::VectorXd{std::move(other)});
}

Eigen::MatrixXd sp::project(const Eigen::MatrixXd & other)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{sp} \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = static_cast<size_t>(2*std::floor(std::min(other.rows(), other.cols())/2));
    const Eigen::MatrixXd square_mat = other(Eigen::seqN(0, shape), Eigen::seqN(0, shape));

    const size_t half_shape = shape / 2;
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(shape, shape);
    J.block(0, half_shape, half_shape, half_shape) = -Eigen::MatrixXd::Identity(half_shape, half_shape);
    J.block(half_shape, 0, half_shape, half_shape) = Eigen::MatrixXd::Identity(half_shape, half_shape);
    Eigen::MatrixXd temp = -J*square_mat;
    return J*(temp + temp.transpose())/2.0;
}

std::ostream & operator<<(std::ostream & os, const sp & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    os << static_cast<const Eigen::MatrixXd>(other.data);
    return os;
}

}

