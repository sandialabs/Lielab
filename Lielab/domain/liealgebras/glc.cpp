#include "glc.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

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

// glc::~glc()

glc::glc(const glc::matrix_t& matrix)
{
    /*! \f{eqnarray*}{(\mathbb{C}^{n \times n}) &\rightarrow& \mathfrak{glc} \f}
    *
    * Constructor instantiating an \f$\mathfrak{glc}\f$ object from an
    * \f$n \times n\f$ imaginary matrix.
    *
    * @param[in] other The object to instantiate from as an imaginary matrix.
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->point.noalias() = matrix;
}

glc glc::basis(const int index, const int shape) 
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{glc} \f}
    *
    * Returns the i'th basis element of the glc algebra.
    * 
    * @param[in] index The basis vector.
    * @param[in] shape The shape of the algebra.
    * @param[out] out The glc element.
    */

    glc out(shape);

    if (index < 0) return out;
    if (index >= out.get_dimension()) return out;

    const int row = (index / 2) / shape;
    const int col = (index / 2) % shape;

    if ((index % 2) == 0)
    {
        out.point(row, col).real(1.0);
    }
    else
    {
        out.point(row, col).imag(1.0);
    }

    return out;
}

glc glc::zero(const int shape)
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

glc glc::from_vector(const Eigen::VectorXd& other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{glc} \f}
    *
    * Constructor instantiating a \f$\mathfrak{glc}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    const int len = static_cast<int>(other.size());
    const int shape = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(len) / 2.0)));
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

glc glc::project(const glc::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in \mathfrak{glc} \f}
    */

    const size_t shape = std::min(matrix.rows(), matrix.cols());
    return glc(matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape)));
}

glc::glc(const int n)
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

    this->point.noalias() = Eigen::MatrixXcd::Zero(n, n);
}

glc glc::from_complex_vector(const Eigen::VectorXcd& other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{glc} \f}
    *
    * Constructor instantiating a \f$\mathfrak{glc}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    const size_t len = other.size();

    const int shape = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(len))));

    Eigen::VectorXd temp = Eigen::VectorXd::Zero(2*shape*shape);

    for (size_t ii = 0; ii < len; ii++)
    {
        temp(2*ii) = other(ii).real();
        temp(2*ii+1) = other(ii).imag();
    }

    glc out(shape);
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

std::string glc::to_string() const
{
    return "gl(" + std::to_string(this->get_shape()) + ", C)";
}

int glc::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    return 2*static_cast<int>(std::pow(this->get_shape(), 2));
}

int glc::get_size() const
{
    return this->get_dimension();
}

bool glc::is_abelian() const
{
    return false;
}

int glc::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the algebra.
    */

    return static_cast<int>(this->point.rows());
}

glc::point_t glc::get_point() const
{
    /*!
    */

    return this->point;
}

Eigen::VectorXd glc::serialize() const
{
    return this->get_vector();
}

void glc::unserialize(const Eigen::VectorXd& serialized)
{
    this->set_vector(serialized);
}

void glc::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
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

    return this->point.cast<std::complex<double>>();
}

Eigen::VectorXd glc::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */
    
    const int dim = this->get_dimension();
    const Eigen::MatrixXcd A = this->get_matrix();

    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);
    int kk = 0;
    for (int ii = 0; ii < this->get_shape(); ii++)
    {
        for (int jj = 0; jj < this->get_shape(); jj++)
        {
            out(kk) = std::real(A(ii, jj));
            out(kk+1) = std::imag(A(ii, jj));
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

    const int vdim = static_cast<int>(vector.size());
    const int max_ind = std::min(this->get_dimension(), vdim);
    
    for (int vind = 0; vind < max_ind; vind++)
    {
        const int rem = vind % 2;
        const int row = (vind / 2) / this->get_shape();
        const int col = (vind / 2) % this->get_shape();
        
        if (rem == 0)
        {
            this->point(row, col).real(vector(vind));
        }
        else
        {
            this->point(row, col).imag(vector(vind));
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

double glc::operator()(const int index) const
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
    
    const int row = (_index/2) / shape;
    const int col = (_index/2) % shape;

    if ((_index % 2) == 0)
    {
        return this->point(row, col).real();
    }

    return this->point(row, col).imag();
}

glc::field_t glc::operator()(const int index1, const int index2) const
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

glc glc::operator+(const glc& other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return glc(this->point + other.point);
}

glc& glc::operator+=(const glc& other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point += other.point;
    return *this;
}

glc glc::operator-(const glc& other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return glc(this->point - other.point);
}

glc& glc::operator-=(const glc& other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point -= other.point;
    return *this;
}

glc glc::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Unary negative of the vector.
    */

    return glc(-this->point);
}

glc glc::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{R}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar product.
    */

    const Eigen::MatrixXcd this_matrix = this->get_matrix();

    return Eigen::MatrixXcd(this_matrix * other);
}

glc operator*(const double other, const glc& rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

glc& glc::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{R}) \rightarrow \mathfrak{glc} \f}
    *
    * In place scalar product.
    */

    this->point *= other;
    return *this;
}

glc glc::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{R}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar division.
    */

    const Eigen::MatrixXcd this_matrix = this->get_matrix();
    return Eigen::MatrixXcd(this_matrix / other);
}

glc& glc::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{R}) \rightarrow \mathfrak{glc} \f}
    *
    * In place scalar division.
    */

    this->point /= other;
    return *this;
}

glc glc::operator*(const std::complex<double> other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar product.
    */

    const Eigen::MatrixXcd this_matrix = this->get_matrix();

    return Eigen::MatrixXcd(this_matrix * other);
}

glc operator*(const std::complex<double> other, const glc& rhs)
{
    /*! \f{equation*}{ (\mathbb{C}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

glc& glc::operator*=(const std::complex<double> other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
    *
    * In place scalar product.
    */

    this->point *= other;
    return *this;
}

glc glc::operator/(const std::complex<double> other) const
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar division.
    */

    const Eigen::MatrixXcd this_matrix = this->get_matrix();
    return Eigen::MatrixXcd(this_matrix / other);
}

glc& glc::operator/=(const std::complex<double> other)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
    *
    * In place scalar division.
    */

    this->point /= other;
    return *this;
}

}
