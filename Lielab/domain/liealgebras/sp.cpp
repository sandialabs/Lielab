#include "sp.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

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

// sp::~sp();

sp::sp(const sp::matrix_t& matrix)
{
    /*! \f{equation}{ (\mathbb{R}^{n \times n}) \rightarrow \mathfrak{sp} \f}
    *
    * Constructor instantiating an \f$\mathfrak{sp}\f$ object from an
    * \f$n \times n\f$ real matrix.
    *
    * @param[in] other The object to instantiate from.
    */

    lielab_assert(matrix.rows() % 2 == 0, "Input matrix shape must be even dimensional.");
    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->point.noalias() = matrix;
}

sp sp::basis(const int index, const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{sp} \f}
    *
    * Returns the i'th basis element of the sp algebra.
    * 
    * @param[in] index The basis vector.
    * @param[in] shape The shape of the algebra.
    * @param[out] out The rn element.
    */

    sp out(shape);
    if (index < 0) return out;

    const int dim = out.get_dimension();
    if (index >= dim) return out;

    Eigen::VectorXd v = Eigen::VectorXd::Zero(dim);
    v(index) = 1.0;
    out.set_vector(v);
    return out;
}

sp sp::zero(const int shape)
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

sp sp::from_vector(const Eigen::VectorXd& other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{sp} \f}
    *
    * Constructor instantiating an \f$\mathfrak{sp}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    const int len = static_cast<int>(other.size());
    int shape = static_cast<int>(std::ceil((-0.5 + std::sqrt(0.25 + 2*len))/(2*0.5)));
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

sp sp::project(const sp::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{sp} \f}
    *
    */

    const int shape = static_cast<int>(2*std::floor(std::min(matrix.rows(), matrix.cols())/2));

    const Eigen::MatrixXd matrix_sq = matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape));

    const int half_shape = shape / 2;
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(shape, shape);
    J.block(0, half_shape, half_shape, half_shape) = -Eigen::MatrixXd::Identity(half_shape, half_shape);
    J.block(half_shape, 0, half_shape, half_shape) = Eigen::MatrixXd::Identity(half_shape, half_shape);
    Eigen::MatrixXd temp = -J*matrix_sq;
    return sp(J*(temp + temp.transpose())/2.0);
}

sp::sp(const int n)
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

    lielab_assert(n % 2 == 0, "Input shape must be even dimensional.");
    this->point.noalias() = Eigen::MatrixXd::Zero(n, n);
}

std::string sp::to_string() const
{
    return "sp(" + std::to_string(this->get_shape()) + ", R)";
}

int sp::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    return this->get_shape() * (this->get_shape() + 1) / 2;
}

int sp::get_size() const
{
    return this->get_dimension();
}

bool sp::is_abelian() const
{
    return false;
}

int sp::get_shape() const
{
    return static_cast<int>(this->point.rows());
}

sp::point_t sp::get_point() const
{
    /*!
    */

    return this->point;
}

Eigen::VectorXd sp::serialize() const
{
    return this->get_vector();
}

void sp::unserialize(const Eigen::VectorXd& serialized)
{
    this->set_vector(serialized);
}

void sp::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
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

    return this->point;
}

Eigen::VectorXd sp::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */

    const int dim = this->get_dimension();
    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);

    if (dim == 0) return out;

    int k = 0;

    for (int ii = 0; ii < this->get_shape()/2; ii++)
    {
        for (int jj = 0; jj < this->get_shape()/2; jj++)
        {
            out(k) = (this->point(ii, jj) - this->point(this->get_shape()/2 + jj, this->get_shape()/2 + ii))/2.0;
            k++;
        }
    }

    for (int ii = 0; ii < this->get_shape()/2; ii++)
    {
        for (int jj = ii; jj < this->get_shape()/2; jj++)
        {
            out(k) = (this->point(ii, this->get_shape()/2 + jj) + this->point(jj, this->get_shape()/2 + ii))/2.0;
            k++;
        }
    }

    for (int ii = 0; ii < this->get_shape()/2; ii++)
    {
        for (int jj = ii; jj < this->get_shape()/2; jj++)
        {
            out(k) = (this->point(this->get_shape()/2 + ii, jj) + this->point(this->get_shape()/2 + jj, ii))/2.0;
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

    const int vdim = static_cast<int>(vector.size());
    const int max_ind = std::min(this->get_dimension(), vdim);

    int k = 0;
    if (k >= max_ind) return;

    for (int ii = 0; ii < this->get_shape()/2; ii++)
    {
        for (int jj = 0; jj < this->get_shape()/2; jj++)
        {
            this->point(ii, jj) = vector(k);
            this->point(this->get_shape()/2 + jj, this->get_shape()/2 + ii) = -vector(k);
            k++;
            if (k >= max_ind) return;
        }
    }

    for (int ii = 0; ii < this->get_shape()/2; ii++)
    {
        for (int jj = ii; jj < this->get_shape()/2; jj++)
        {
            this->point(ii, this->get_shape()/2 + jj) = vector(k);
            this->point(jj, this->get_shape()/2 + ii) = vector(k);
            k++;
            if (k >= max_ind) return;
        }
    }

    for (int ii = 0; ii < this->get_shape()/2; ii++)
    {
        for (int jj = ii; jj < this->get_shape()/2; jj++)
        {
            this->point(this->get_shape()/2 + ii, jj) = vector(k);
            this->point(this->get_shape()/2 + jj, ii) = vector(k);
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

double sp::operator()(const int index) const
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

    const Eigen::VectorXd vector = this->get_vector();
    return vector(_index);
}

sp::field_t sp::operator()(const int index1, const int index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int shape = this->get_shape();
    if (shape == 0) return nan;
    if (shape % 2 == 1) return nan;

    // If input index is negative, index from the back of the array
    const int _index1 = (index1 < 0) ? shape + index1 : index1;
    const int _index2 = (index2 < 0) ? shape + index2 : index2;

    // Error check for out of bounds
    if (_index1 < 0) return nan;
    if (_index1 >= shape) return nan;
    if (_index2 < 0) return nan;
    if (_index2 >= shape) return nan;

    return this->point(_index1, _index2);
}

sp sp::operator+(const sp& other) const
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * Addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return sp(this->point + other.point);
}

sp& sp::operator+=(const sp& other)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point += other.point;
    return *this;
}

sp sp::operator-(const sp& other) const
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return sp(this->point - other.point);
}

sp& sp::operator-=(const sp& other)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point -= other.point;
    return *this;
}

sp sp::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * Unary negative of the vector.
    */

    return sp(-this->point);
}

sp sp::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
    *
    * Scalar product.
    */

    Eigen::MatrixXd out = this->point * other;
    return out;
}

sp operator*(const double other, const sp& rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

sp& sp::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
    *
    * In place scalar product.
    */

    this->point *= other;
    return *this;
}

sp sp::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
    *
    * Scalar division.
    */

    Eigen::MatrixXd out = this->point / other;
    return out;
}

sp& sp::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
    *
    * In place scalar division.
    */

    this->point /= other;
    return *this;
}

}

