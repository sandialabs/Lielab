#include "su.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string su::to_string() const
{
    return "su(" + std::to_string(this->_shape) + ")";
}

su::su() : su(0)
{
    /*! \f{equation*}{ () \rightarrow \mathfrak{su} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::su x, y, z;
    * 
    */

}

su::su(const size_t n)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{su} \f}
    *
    * Constructor instantiating an \f$\mathfrak{su}\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::su x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */
    
    this->_shape = n;
    this->data = Eigen::MatrixXcd::Zero(n, n);
}

su su::basis(const ptrdiff_t i, const size_t n)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{su} \f}
    *
    * Returns the i'th basis element of the su algebra.
    * 
    * @param[in] i The basis vector.
    * @param[in] n The size of the algebra.
    * @param[out] out The su element.
    */

    su out(n);
    if (i < 0) return out;

    const size_t ind = static_cast<size_t>(i);

    const size_t dim = out.get_dimension();
    if (ind >= dim) return out;

    Eigen::VectorXd v = Eigen::VectorXd::Zero(dim);
    v(ind) = 1.0;
    out.set_vector(v);
    return out;
}

su su::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{su} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The su element. 
    */

    return su(shape);
}

size_t su::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    if (this->_shape == 0) return 0;

    return this->_shape * this->_shape - 1;
}

Eigen::VectorXd su::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    *
    * Sources:
    *     - Georgi, Howard. Lie algebras in particle physics: from isospin to unified theories.
    *         Taylor & Francis, 2000.
    *     - Stover, Christopher. "Generalized Gell-Mann Matrix." From MathWorld--A Wolfram Web Resource,
    *         created by Eric W. Weisstein. https://mathworld.wolfram.com/GeneralizedGell-MannMatrix.html 
    */

    if (this->_shape <= 1)
    {
        return Eigen::VectorXd::Zero(0);
    }

    const size_t dim = this->get_dimension();
    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);

    // General case of su(2+). Use's the Generalized Gell-Mann matrices
    size_t k = 0;

    // Symmetric
    for (size_t jj = 1; jj < this->_shape; jj++)
    {
        for (size_t ii = 0; ii < jj; ii++)
        {
            out(k) = std::imag(this->data(ii, jj) + this->data(jj, ii))/2.0;
            k++;
        }
    }

    // Anti-symmetric
    for (size_t jj = 1; jj < this->_shape; jj++)
    {
        for (size_t ii = 0; ii < jj; ii++)
        {
            out(k) = std::real(this->data(jj, ii) - this->data(ii, jj))/2.0;
            k++;
        }
    }

    // Diagonal
    size_t zz = this->_shape;
    k = out.size() - 1;
    Eigen::MatrixXcd temp = this->get_matrix();
    for (size_t yy = this->_shape - 1; yy >= 1; yy--)
    {
        const double multiplier = std::sqrt(2.0/((zz-1)*(zz)));
        out(k) = -std::imag(temp(yy, yy))/(multiplier*(zz - 1));

        for (size_t ii = 0; ii < yy; ii++)
        {
            temp(ii, ii) -= std::complex<double>(0.0, multiplier*out(k)); // Do not use std::imag(). This doesn't work w/ inplace operations with Eigen.
        }

        zz--;
        k--;
    }

    return out;
}

void su::set_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ \mathfrak{su} := \mathbb{R}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    * 
    * Sources:
    *     - Georgi, Howard. Lie algebras in particle physics: from isospin to unified theories.
    *         Taylor & Francis, 2000.
    *     - Stover, Christopher. "Generalized Gell-Mann Matrix." From MathWorld--A Wolfram Web Resource,
    *         created by Eric W. Weisstein. https://mathworld.wolfram.com/GeneralizedGell-MannMatrix.html 
    */
    
    const size_t vdim = vector.size();
    const size_t dim = this->get_dimension();

    // su(0) and su(1) are 0-dimensional. Do nothing.
    if (this->_shape <= 1) return;

    // Reuse data from current object if vdim < dim
    // TODO: This function could be made more efficient if
    //       it "removed" the data from each index, then assigned the
    //       new values instead of rewriting the entire matrix.
    const Eigen::VectorXd vec0 = this->get_vector();
    Eigen::VectorXd vec_assign = Eigen::VectorXd::Zero(dim);
    for (size_t ii = 0; ii < dim; ii++)
    {
        if (ii < vdim)
        {
            vec_assign(ii) = vector(ii);
        }
        else
        {
            vec_assign(ii) = vec0(ii);
        }
    }
    this->data = data_t::Zero(this->_shape, this->_shape);

    // General case of su(2+). Use's the Generalized Gell-Mann matrices
    size_t k = 0;

    // Symmetric
    for (size_t jj = 1; jj < this->_shape; jj++)
    {
        for (size_t ii = 0; ii < jj; ii++)
        {
            this->data(ii, jj) = std::complex<double>(0.0, vec_assign(k));
            this->data(jj, ii) = std::complex<double>(0.0, vec_assign(k));
            k++;
        }
    }

    // Anti-symmetric
    for (size_t jj = 1; jj < this->_shape; jj++)
    {
        for (size_t ii = 0; ii < jj; ii++)
        {
            this->data(ii, jj) += std::complex<double>(-vec_assign(k), 0.0);
            this->data(jj, ii) += std::complex<double>(vec_assign(k), 0.0);
            k++;
        }
    }

    // Diagonal
    size_t zz = 2;
    while (k < dim)
    {
        const double multiplier = std::sqrt(2.0/((zz-1)*(zz)));

        for (size_t ii = 0; ii < zz; ii++)
        {
            if (ii == (zz - 1))
            {
                this->data(ii, ii) -= std::complex<double>(0.0, multiplier*(zz - 1)*vec_assign(k));
            }
            else
            {
                this->data(ii, ii) += std::complex<double>(0.0, multiplier*vec_assign(k));
            }
        }

        zz++;
        k++;
    }
}

void su::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

su::matrix_t su::get_matrix() const
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

double su::operator()(const ptrdiff_t index) const
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

    return this->get_vector()(_index);
}

std::complex<double> su::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
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

su su::operator+(const su & other) const
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = this->get_matrix();
    const Eigen::MatrixXcd rhs_matrix = other.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return su(new_matrix);
}

su & su::operator+=(const su & other)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    assert(this->_shape == other.get_shape());
    this->data += other.data;
    return *this;
}

su su::operator-(const su & other) const
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(this->_shape, other.get_shape());
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = this->get_matrix();
    const Eigen::MatrixXcd rhs_matrix = other.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return su(new_matrix);
}

su & su::operator-=(const su & other)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    assert(this->_shape == other.get_shape());
    this->data -= other.data;
    return *this;
}

su su::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * Unary negative of the vector.
    */

    Eigen::MatrixXcd out = -this->data;
    return out;
}

su su::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    *
    * Scalar product.
    */

    Eigen::MatrixXcd out = this->data * other;
    return out;
}

su operator*(const double other, const su & rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

su su::operator*(const std::complex<double> other) const
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{C}) \rightarrow \mathfrak{su} \f}
    *
    * Scalar product.
    */

    Eigen::MatrixXcd out = this->data * other;
    return out;
}

su operator*(const std::complex<double> other, const su & rhs)
{
    /*! \f{equation*}{ (\mathbb{C}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

su & su::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    *
    * In place scalar product.
    */

    this->data *= other;
    return *this;
}

su & su::operator*=(const std::complex<double> other)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{C}) \rightarrow \mathfrak{su} \f}
    *
    * In place scalar product.
    */

    this->data *= other;
    return *this;
}

su su::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    *
    * Scalar division.
    */

    Eigen::MatrixXcd out = this->data / other;
    return out;
}

su su::operator/(const std::complex<double> other) const
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{C}) \rightarrow \mathfrak{su} \f}
    *
    * Scalar division.
    */

    Eigen::MatrixXcd out = this->data / other;
    return out;
}

su & su::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    *
    * In place scalar division.
    */

    this->data /= other;
    return *this;
}

su & su::operator/=(const std::complex<double> other)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{C}) \rightarrow \mathfrak{su} \f}
    *
    * In place scalar division.
    */

    this->data /= other;
    return *this;
}

su su::from_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{su} \f}
    *
    * Constructor instantiating an \f$\mathfrak{su}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    const size_t len = vector.size();
    const size_t shape = static_cast<size_t>(std::ceil((std::sqrt(4*(len+1)))/2));
    
    su out(shape);
    out.set_vector(vector);
    return out;
}

su su::from_vector(std::initializer_list<double> vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{su} \f}
    *
    * Constructor instantiating an \f$\mathfrak{su}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] vector The object to instantiate from as a real vector.
    */

    return su::from_vector(Eigen::VectorXd{std::move(vector)});
}

// TODO: Project function

std::ostream & operator<<(std::ostream & os, const su & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    os << static_cast<const Eigen::MatrixXcd>(other.data);
    return os;
}

}
