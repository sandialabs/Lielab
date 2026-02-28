#include "cn.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <complex>

namespace Lielab::domain
{

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

// cn::~cn()

cn::cn(const cn::matrix_t& matrix)
{
    /*! \f{eqnarray*}{(\mathbb{C}^{n \times n}) &\rightarrow& \mathfrak{cn} \\ (\mathbb{C}^{n \times 1}) &\rightarrow& \mathfrak{cn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{cn}\f$ object from either an
    * \f$n \times n\f$ imaginary matrix or \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary matrix.
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->_shape = static_cast<int>(matrix.rows());

    if (this->_shape == 0)
    {
        this->point = Eigen::VectorXcd::Zero(0);
        return;
    }

    this->point = Eigen::VectorXcd::Zero(this->_shape - 1);
    for (int ii = 0; ii < this->_shape - 1; ii++)
    {
        this->point(ii) = matrix(ii, this->_shape - 1);
    }
}

cn cn::basis(const int index, const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{cn} \f}
    *
    * Returns the i'th basis element of the cn algebra.
    * 
    * @param[in] index The basis vector.
    * @param[in] shape The shape of the algebra.
    * @param[out] out The cn element.
    */

    cn out = cn::zero(shape);
    if (index < 0) return out;

    const int dim = out.get_dimension();

    if (index >= dim) return out;

    const size_t indz = index / 2;
    const size_t rem = index % 2;
    if (rem == 0) out.point(indz) = std::complex<double>(1.0, 0.0);
    if (rem == 1) out.point(indz) = std::complex<double>(0.0, 1.0);

    return out;
}

cn cn::zero(const int shape)
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

cn cn::from_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{cn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{cn}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] vector The object to instantiate from as an imaginary vector.
    */

    const int shape = static_cast<int>(std::ceil(vector.size()/2.0)) + 1;
    cn out = cn::zero(shape);
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

cn cn::project(const cn::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in \mathfrak{cn} \f}
    *
    */

    const size_t shape = std::min(matrix.rows(), matrix.cols());
    return cn(matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape)));
}

cn::cn(const int n)
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
    this->point.noalias() = Eigen::VectorXcd::Zero(n);
}

cn cn::from_complex_vector(const Eigen::VectorXcd& other)
{
    /*! \f{equation}{(\mathbb{C}^{n \times 1}) \rightarrow \mathfrak{cn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{cn}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    const int n = static_cast<int>(other.size());
    cn out(n);
    out.point = other;
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

std::string cn::to_string() const
{
    const size_t shape = this->get_shape();
    if (shape == 0) return "c^nan";
    return "c^" + std::to_string(shape-1);
}

int cn::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    return 2*static_cast<int>(this->point.size());
}

int cn::get_size() const
{
    return this->get_dimension();
}

bool cn::is_abelian() const
{
    return true;
}

int cn::get_shape() const
{
    return this->_shape;
}

cn::point_t cn::get_point() const
{
    /*!
    */

    return this->point;
}

Eigen::VectorXd cn::serialize() const
{
    return this->get_vector();
}

void cn::unserialize(const Eigen::VectorXd& serialized)
{
    this->set_vector(serialized);
}

void cn::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
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

    for (int ii = 0; ii < this->_shape - 1; ii++)
    {
        out(ii, this->_shape - 1) = this->point(ii);
    }

    return out;
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
        out(kk) = std::real(this->point(ii));
        kk += 1;
        out(kk) = std::imag(this->point(ii));
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

    const int vdim = static_cast<int>(vector.size());
    const int max_ind = std::min(this->get_dimension(), vdim);

    int dind = 0;
    for (int vind = 0; vind < max_ind; vind++)
    {
        const int rem = vind % 2;

        if (rem == 0)
        {
            this->point(dind).real(vector(vind));
        }
        else
        {
            this->point(dind).imag(vector(vind));
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

double cn::operator()(const int index) const
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

    if (_index % 2 == 0) return this->point(_index/2).real();
    return this->point(_index/2).imag();
}

cn::field_t cn::operator()(const int index1, const int index2) const
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

    if (_index1 == shape - 1) return std::complex<double>(0.0, 0.0);
    if (_index2 != shape - 1) return std::complex<double>(0.0, 0.0);

    return this->point(_index1);
}

Eigen::VectorXcd cn::to_complex_vector() const
{
    /*
     * Returns complex vector.
     */
    
    return this->point;
}

const cn::field_t& cn::operator[](const int index) const
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    */

    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim/2 + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < dim/2), "Index " + std::to_string(index) + " is out of bounds for cn of length " + std::to_string(dim/2));

    return this->point(_index);
}

cn::field_t& cn::operator[](const int index)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    */

    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim/2 + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < dim/2), "Index " + std::to_string(index) + " is out of bounds for cn of length " + std::to_string(dim/2));

    return this->point(_index);
}

cn cn::operator+(const cn& other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return cn::from_complex_vector(this->point + other.point);
}

cn& cn::operator+=(const cn& other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point += other.point;
    return *this;
}

cn cn::operator-(const cn& other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return cn::from_complex_vector(this->point - other.point);
}

cn& cn::operator-=(const cn& other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point -= other.point;
    return *this;
}

cn cn::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Unary negative of the vector.
    */

    return cn::from_complex_vector(-this->point);
}

cn cn::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */

    return cn::from_complex_vector(this->point*other);
}

cn operator*(const double other, const cn& rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */

    return cn::from_complex_vector(other*rhs.point);
}

cn& cn::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar product.
    */

    this->point *= other;
    return *this;
}

cn cn::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar division.
    */

    return cn::from_complex_vector(this->point/other);
}

cn& cn::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar division.
    */

    this->point /= other;
    return *this;
}

cn cn::operator*(const std::complex<int> other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */
    
    const std::complex<double> otherd = std::complex<double>(static_cast<double>(other.real()),
                                                             static_cast<double>(other.imag()));
    return cn::from_complex_vector(this->point*otherd);
}

cn cn::operator*(const std::complex<double> other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */

    return cn::from_complex_vector(this->point*other);
}

cn operator*(const std::complex<int> other, const cn& rhs)
{
    /*! \f{equation*}{ (\mathbb{C}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */
    
    const std::complex<double> otherd = std::complex<double>(static_cast<double>(other.real()),
                                                             static_cast<double>(other.imag()));
    return cn::from_complex_vector(otherd*rhs.point);
}

cn operator*(const std::complex<double> other, const cn& rhs)
{
    /*! \f{equation*}{ (\mathbb{C}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar product.
    */

    return cn::from_complex_vector(other*rhs.point);
}

cn& cn::operator*=(const std::complex<int> other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar product.
    */
    
    const std::complex<double> otherd = std::complex<double>(static_cast<double>(other.real()),
                                                             static_cast<double>(other.imag()));
    this->point *= otherd;
    return *this;
}

cn& cn::operator*=(const std::complex<double> other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar product.
    */

    this->point *= other;
    return *this;
}

cn cn::operator/(const std::complex<int> other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar division.
    */

    const std::complex<double> otherd = std::complex<double>(static_cast<double>(other.real()),
                                                             static_cast<double>(other.imag()));
    return cn::from_complex_vector(this->point/otherd);
}

cn cn::operator/(const std::complex<double> other) const
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * Scalar division.
    */

    return cn::from_complex_vector(this->point/other);
}

cn& cn::operator/=(const std::complex<int> other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar division.
    */

    const std::complex<double> otherd = std::complex<double>(static_cast<double>(other.real()),
                                                             static_cast<double>(other.imag()));
    this->point /= otherd;
    return *this;
}

cn& cn::operator/=(const std::complex<double> other)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
    *
    * In place scalar division.
    */

    this->point /= other;
    return *this;
}

}
