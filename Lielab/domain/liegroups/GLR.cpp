#include "GLR.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>


namespace Lielab::domain
{

GLR::GLR() : GLR(0)
{
    /*! \f{equation*}{ () \rightarrow GLR \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::GLR x, y, z;
    * 
    */

}

// GLR::~GLR();

GLR::GLR(const GLR::matrix_t& matrix)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& GLR \\ (\mathbb{R}^{n \times 1}) &\rightarrow& GLR \f}
    *
    * Constructor instantiating an \f$GLR\f$ object from either an
    * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->point.noalias() = matrix;
}

GLR GLR::identity(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{GLR} \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The GLR element.
    */

    return GLR(shape);
}

GLR GLR::project(const GLR::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in GLR \f}
    *
    */

    const size_t shape = std::min(matrix.rows(), matrix.cols());
    return GLR(matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape)));
}

GLR::GLR(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow GLR \f}
    *
    * Constructor instantiating an \f$GLR\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::GLR x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */
    
    this->point.noalias() = point_t::Identity(shape, shape);
}

std::string GLR::to_string() const
{
    return "GL(" + std::to_string(this->get_shape()) + ", R)";
}

int GLR::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    return static_cast<int>(std::pow(this->get_shape(), 2));
}

int GLR::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
    *
    * Gets the size of the data representation.
    */

    return static_cast<int>(std::pow(this->get_shape(), 2));
}

bool GLR::is_abelian() const
{
    return false;
}

int GLR::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return static_cast<int>(this->point.rows());
}

GLR::point_t GLR::get_point() const
{
    return this->point;
}

Eigen::VectorXd GLR::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */

    return this->point.reshaped<Eigen::RowMajor>();
}

void GLR::unserialize(const Eigen::VectorXd& serialized)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the GLR object from a serialized vector.
    */

    const int vdim = static_cast<int>(serialized.size());
    const int max_ind = std::min(this->get_dimension(), vdim);

    for (int vind = 0; vind < max_ind; vind++)
    {
        const int row = vind / this->get_shape();
        const int col = vind % this->get_shape();
        this->point(row, col) = serialized(vind);
    }
}

void GLR::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

GLR::matrix_t GLR::get_matrix() const
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

GLR::field_t GLR::operator()(const int index1, const int index2) const
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

GLR GLR::operator*(const GLR& other) const
{
    /*! \f{equation*}{ (GLR, GLR) \rightarrow GLR \f}
    *
    * Group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return GLR(this->point * other.point);
}

GLR& GLR::operator*=(const GLR& other)
{
    /*! \f{equation*}{ (GLR, GLR) \rightarrow GLR \f}
    *
    * In place group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point *= other.point;
    return *this;
}

GLR GLR::inverse() const
{
    /*! \f{equation*}{ (GLR) \rightarrow GLR \f}
    * 
    * Returns the inverse.
    */

    return GLR(this->point.inverse());
}

}
