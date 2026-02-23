#include "su.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

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

// su::~su();

su::su(const su::matrix_t& matrix)
{
    /*! \f{equation}{ (\mathbb{C}^{n \times n}) \rightarrow \mathfrak{su} \f}
    *
    * Constructor instantiating an \f$\mathfrak{su}\f$ object from an
    * \f$n \times n\f$ real matrix.
    */
    
    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->point.noalias() = matrix;
}

su su::basis(const int index, const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{su} \f}
    *
    * Returns the i'th basis element of the su algebra.
    * 
    * @param[in] index The basis vector.
    * @param[in] shape The shape of the algebra.
    * @param[out] out The su element.
    */

    su out(shape);
    if (index < 0) return out;

    const int dim = out.get_dimension();
    if (index >= dim) return out;

    Eigen::VectorXd v = Eigen::VectorXd::Zero(dim);
    v(index) = 1.0;
    out.set_vector(v);
    return out;
}

su su::zero(const int shape)
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

su su::from_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{su} \f}
    *
    * Constructor instantiating an \f$\mathfrak{su}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    const int len = static_cast<int>(vector.size());
    const int shape = static_cast<int>(std::ceil((std::sqrt(4*(len+1)))/2));
    
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

su su::project(const su::matrix_t& matrix)
{
    const size_t shape = std::min(matrix.rows(), matrix.cols());

    if (shape == 0) return su::zero(0);

    su::matrix_t matrix_sq = matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape));

    const su::matrix_t temp = (matrix_sq - matrix_sq.adjoint())/(std::complex<double>(2.0, 0.0));
    const std::complex<double> trace = temp.trace();
    const su::matrix_t matrix_projected = temp - trace/static_cast<double>(shape)*Eigen::MatrixXcd::Identity(shape, shape);
    return su(matrix_projected);
}

su::su(const int n)
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
    
    this->point.noalias() = Eigen::MatrixXcd::Zero(n, n);
}

std::string su::to_string() const
{
    return "su(" + std::to_string(this->get_shape()) + ")";
}

int su::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    if (this->get_shape() == 0) return 0;

    return this->get_shape() * this->get_shape() - 1;
}

int su::get_size() const
{
    return this->get_dimension();
}

bool su::is_abelian() const
{
    return false;
}

int su::get_shape() const
{
    return static_cast<int>(this->point.rows());
}

su::point_t su::get_point() const
{
    /*!
    */

    return this->point;
}

Eigen::VectorXd su::serialize() const
{
    return this->get_vector();
}

void su::unserialize(const Eigen::VectorXd& serialized)
{
    this->set_vector(serialized);
}

void su::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
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

    return this->point;
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

    if (this->get_shape() <= 1)
    {
        return Eigen::VectorXd::Zero(0);
    }

    const size_t dim = this->get_dimension();
    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);

    // General case of su(2+). Use's the Generalized Gell-Mann matrices
    int k = 0;

    // Symmetric
    for (int jj = 1; jj < this->get_shape(); jj++)
    {
        for (int ii = 0; ii < jj; ii++)
        {
            out(k) = std::imag(this->point(ii, jj) + this->point(jj, ii))/2.0;
            k++;
        }
    }

    // Anti-symmetric
    for (int jj = 1; jj < this->get_shape(); jj++)
    {
        for (int ii = 0; ii < jj; ii++)
        {
            out(k) = std::real(this->point(jj, ii) - this->point(ii, jj))/2.0;
            k++;
        }
    }

    // Diagonal
    int zz = this->get_shape();
    k = static_cast<int>(out.size()) - 1;
    Eigen::MatrixXcd temp = this->get_matrix();
    for (int yy = this->get_shape() - 1; yy >= 1; yy--)
    {
        const double multiplier = std::sqrt(2.0/((zz-1)*(zz)));
        out(k) = -std::imag(temp(yy, yy))/(multiplier*(zz - 1));

        for (int ii = 0; ii < yy; ii++)
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
    
    const int vdim = static_cast<int>(vector.size());
    const int dim = this->get_dimension();

    // su(0) and su(1) are 0-dimensional. Do nothing.
    if (this->get_shape() <= 1) return;

    // Reuse data from current object if vdim < dim
    // TODO: This function could be made more efficient if
    //       it "removed" the data from each index, then assigned the
    //       new values instead of rewriting the entire matrix.
    const Eigen::VectorXd vec0 = this->get_vector();
    Eigen::VectorXd vec_assign = Eigen::VectorXd::Zero(dim);
    for (int ii = 0; ii < dim; ii++)
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
    this->point = point_t::Zero(this->get_shape(), this->get_shape());

    // General case of su(2+). Use's the Generalized Gell-Mann matrices
    int k = 0;

    // Symmetric
    for (int jj = 1; jj < this->get_shape(); jj++)
    {
        for (int ii = 0; ii < jj; ii++)
        {
            this->point(ii, jj) = std::complex<double>(0.0, vec_assign(k));
            this->point(jj, ii) = std::complex<double>(0.0, vec_assign(k));
            k++;
        }
    }

    // Anti-symmetric
    for (int jj = 1; jj < this->get_shape(); jj++)
    {
        for (int ii = 0; ii < jj; ii++)
        {
            this->point(ii, jj) += std::complex<double>(-vec_assign(k), 0.0);
            this->point(jj, ii) += std::complex<double>(vec_assign(k), 0.0);
            k++;
        }
    }

    // Diagonal
    int zz = 2;
    while (k < dim)
    {
        const double multiplier = std::sqrt(2.0/((zz-1)*(zz)));

        for (int ii = 0; ii < zz; ii++)
        {
            if (ii == (zz - 1))
            {
                this->point(ii, ii) -= std::complex<double>(0.0, multiplier*(zz - 1)*vec_assign(k));
            }
            else
            {
                this->point(ii, ii) += std::complex<double>(0.0, multiplier*vec_assign(k));
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

double su::operator()(const int index) const
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the column vector representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim + index : index;

    // Error check for out of bounds
    if (_index < 0) return nan;
    if (_index >= dim) return nan;

    return this->get_vector()(_index);
}

su::field_t su::operator()(const int index1, const int index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int shape = this->get_shape();
    if (shape == 0) return std::complex<double>(nan, nan);

    // If input index is negative, index from the back of the array
    const int _index1 = (index1 < 0) ? shape + index1 : index1;
    const int _index2 = (index2 < 0) ? shape + index2 : index2;

    // Error check for out of bounds
    if (_index1 < 0) return std::complex<double>(nan, nan);
    if (_index1 >= shape) return std::complex<double>(nan, nan);
    if (_index2 < 0) return std::complex<double>(nan, nan);
    if (_index2 >= shape) return std::complex<double>(nan, nan);

    return this->point(_index1, _index2);
}

su su::operator+(const su& other) const
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * Addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return su(this->point + other.point);
}

su& su::operator+=(const su& other)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point += other.point;
    return *this;
}

su su::operator-(const su& other) const
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return su(this->point - other.point);
}

su& su::operator-=(const su& other)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point -= other.point;
    return *this;
}

su su::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * Unary negative of the vector.
    */

    return su(-this->point);
}

su su::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    *
    * Scalar product.
    */

    Eigen::MatrixXcd out = this->point * other;
    return out;
}

su operator*(const double other, const su& rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

su& su::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    *
    * In place scalar product.
    */

    this->point *= other;
    return *this;
}

su su::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    *
    * Scalar division.
    */

    Eigen::MatrixXcd out = this->point / other;
    return out;
}

su& su::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    *
    * In place scalar division.
    */

    this->point /= other;
    return *this;
}

}
