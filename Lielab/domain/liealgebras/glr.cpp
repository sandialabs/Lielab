#include "glr.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

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

// glr::~glr();

glr::glr(const glr::matrix_t& matrix)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{glr} \\ (\mathbb{R}^{n \times 1}) &\rightarrow& \mathfrak{glr} \f}
    *
    * Constructor instantiating an \f$\mathfrak{glr}\f$ object from either an
    * \f$n \times n\f$ real matrix.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->point.noalias() = matrix;
}

glr glr::basis(const int index, const int shape) 
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{glr} \f}
    *
    * Returns the i'th basis element of the glr algebra.
    * 
    * @param[in] index The basis vector.
    * @param[in] shape The shape of the algebra.
    * @param[out] out The glr element.
    */

    glr out(shape);

    if (index < 0) return out;
    if (index >= out.get_dimension()) return out;

    const int row = index / shape;
    const int col = index % shape;

    out.point(row, col) = 1.0;

    return out;
}

glr glr::zero(const int shape)
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

glr glr::from_vector(const Eigen::VectorXd& other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{glr} \f}
    *
    * Constructor instantiating a \f$\mathfrak{glr}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    const int len = static_cast<int>(other.size());
    const int shape = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(len))));
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

glr glr::project(const glr::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{glr} \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(matrix.rows(), matrix.cols());
    return glr(matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape)));
}

glr::glr(const int n)
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

    this->point.noalias() = Eigen::MatrixXd::Zero(n, n);
}

std::string glr::to_string() const
{
    return "gl(" + std::to_string(this->get_shape()) + ", R)";
}

int glr::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    return static_cast<int>(std::pow(this->get_shape(), 2));
}

int glr::get_size() const
{
    return this->get_dimension();
}

bool glr::is_abelian() const
{
    return false;
}

int glr::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the algebra.
    */

    return static_cast<int>(this->point.rows());
}

glr::point_t glr::get_point() const
{
    /*!
    */

    return this->point;
}

Eigen::VectorXd glr::serialize() const
{
    return this->get_vector();
}

void glr::unserialize(const Eigen::VectorXd& serialized)
{
    this->set_vector(serialized);
}

void glr::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
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

    return this->point;
}

Eigen::VectorXd glr::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */
    
    return this->point.reshaped<Eigen::RowMajor>();
}

void glr::set_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ \mathfrak{glr} := \mathbb{R}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    */

    const int vdim = static_cast<int>(vector.size());
    const int max_ind = std::min(this->get_dimension(), vdim);

    for (int vind = 0; vind < max_ind; vind++)
    {
        const int row = vind / this->get_shape();
        const int col = vind % this->get_shape();
        this->point(row, col) = vector(vind);
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

double glr::operator()(const int index) const
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the vector representation.
    */
    
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim + index : index;

    // Error check for out of bounds
    if (_index < 0) return nan;
    if (_index >= dim) return nan;

    const int shape = this->get_shape();

    const int row = _index / shape;
    const int col = _index % shape;

    return this->point(row, col);
}

glr::field_t glr::operator()(const int index1, const int index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */
    
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int shape = this->get_shape();
    if (shape == 0) return nan;

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

glr glr::operator+(const glr& other) const
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * Addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return glr(this->point + other.point);
}

glr& glr::operator+=(const glr& other)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point += other.point;
    return *this;
}

glr glr::operator-(const glr& other) const
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return glr(this->point - other.point);
}

glr& glr::operator-=(const glr& other)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point -= other.point;
    return *this;
}

glr glr::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * Unary negative of the vector.
    */

    return glr(-this->point);
}

glr glr::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    *
    * Scalar product.
    */

    return Eigen::MatrixXd(this->point * other);
}

glr operator*(const double other, const glr& rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{glr}) \rightarrow \mathfrak{glr} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

glr& glr::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    *
    * In place scalar product.
    */

    this->point *= other;
    return *this;
}

glr glr::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    *
    * Scalar division.
    */

    return Eigen::MatrixXd(this->point / other);
}

glr& glr::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    *
    * In place scalar division.
    */

    this->point /= other;
    return *this;
}

}
