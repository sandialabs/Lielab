#include "so.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

so::so() : so(0)
{
    /*! \f{equation*}{() \rightarrow \mathfrak{so} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::so x, y, z;
    * 
    */

}

// so::~so();

so::so(const so::matrix_t& matrix)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{so} \\ (\mathbb{R}^{n \times 1}) &\rightarrow& \mathfrak{so} \f}
    *
    * Constructor instantiating an \f$\mathfrak{so}\f$ object from either an
    * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");
    
    this->point.noalias() = matrix;
}

so so::basis(const int index, const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{so} \f}
    *
    * Returns the i'th basis element of the so algebra.
    * 
    * @param[in] index The basis vector.
    * @param[in] shape The size of the algebra.
    * @param[out] out The so element.
    */

    so out(shape);
    if (index < 0) return out;

    const int dim = out.get_dimension();
    if (index >= dim) return out;

    Eigen::VectorXd v = Eigen::VectorXd::Zero(dim);
    v(index) = 1.0;
    out.set_vector(v);
    return out;
}

so so::zero(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{so} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The so element. 
    */

    return so(shape);
}

so so::from_vector(const Eigen::VectorXd& other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{so} \f}
    *
    * Constructor instantiating an \f$\mathfrak{so}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    const int len = static_cast<int>(other.size());
    const int shape = static_cast<int>(std::ceil(std::sqrt(2.0*len + 0.25) + 0.5));

    so out = so::zero(shape);
    out.set_vector(other);
    return out;
}

so so::from_vector(std::initializer_list<double> other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{so} \f}
    *
    * Constructor instantiating an \f$\mathfrak{so}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    return so::from_vector(Eigen::VectorXd{std::move(other)});
}

so so::project(const so::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{so} \f}
    *
    * Projects a matrix suitable for data.
    */
    
    const size_t shape = std::min(matrix.rows(), matrix.cols());

    const Eigen::MatrixXd matrix_sq = matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape));

    return so((matrix_sq - matrix_sq.transpose())/2.0);
}

so::so(const int n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{so} \f}
    *
    * Constructor instantiating an \f$\mathfrak{so}\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::so x(3), y(4), z(5);
    * 
    * @param[in] shape The shape of the data matrix.
    */
    
    this->point.noalias() = Eigen::MatrixXd::Zero(n, n);
}

std::string so::to_string() const
{
    return "so(" + std::to_string(this->get_shape()) + ")";
}

int so::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    if (this->get_shape() == 0) return 0; // TODO: Return nan?

    return this->get_shape() * (this->get_shape() - 1) / 2;
}

int so::get_size() const
{
    return this->get_dimension();
}

bool so::is_abelian() const
{
    return false;
}

int so::get_shape() const
{
    return static_cast<int>(this->point.rows());
}

so::point_t so::get_point() const
{
    /*!
    */

    return this->point;
}

Eigen::VectorXd so::serialize() const
{
    return this->get_vector();
}

void so::unserialize(const Eigen::VectorXd& serialized)
{
    this->set_vector(serialized);
}

void so::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

so::matrix_t so::get_matrix() const
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

Eigen::VectorXd so::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */

    const size_t dim = this->get_dimension();
    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);
    if (this->get_shape() == 0) return out;

    int k = 0;

    for (size_t ii = this->get_shape() - 1; ii > 0; ii--)
    {
        for (size_t jj = this->get_shape(); jj > ii; jj--)
        {
            out(k) = this->point(ii-1, jj-1)/std::pow(-1.0, ii+jj);
            k++;
        }
    }

    return out;
}

void so::set_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ \mathfrak{so} := \mathbb{R}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    */

    const int vdim = static_cast<int>(vector.size());
    const int max_ind = std::min(this->get_dimension(), vdim);
    int k = 0;
    if (k >= max_ind) return;

    for (int ii = this->get_shape() - 1; ii > 0; ii--)
    {
        for (int jj = this->get_shape(); jj > ii; jj--)
        {
            this->point(ii-1, jj-1) =  std::pow(-1, (ii+jj))*vector(k);
            this->point(jj-1, ii-1) = -std::pow(-1, (ii+jj))*vector(k);
            k++;
            if (k >= max_ind) return;
        }
    }
}

void so::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

double so::operator()(const int index) const
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

so::field_t so::operator()(const int index1, const int index2) const
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

so so::operator+(const so& other) const
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}) \rightarrow \mathfrak{so} \f}
    *
    * Addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return so(this->point + other.point);
}

so& so::operator+=(const so& other)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}) \rightarrow \mathfrak{so} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point += other.point;
    return *this;
}

so so::operator-(const so& other) const
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}) \rightarrow \mathfrak{so} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return so(this->point - other.point);
}

so& so::operator-=(const so& other)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}) \rightarrow \mathfrak{so} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point -= other.point;
    return *this;
}

so so::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{so}) \rightarrow \mathfrak{so} \f}
    *
    * Unary negative of the vector.
    */
    
    return so(-this->point);
}

so so::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{so} \f}
    *
    * Scalar product.
    */

    return Eigen::MatrixXd(this->point * other);
}

so operator*(const double other, const so& rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{so}) \rightarrow \mathfrak{so} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

so& so::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{so} \f}
    *
    * In place scalar multiplication.
    */

    this->point *= other;
    return *this;
}

so so::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{so} \f}
    *
    * Scalar division.
    */

    Eigen::MatrixXd out = this->point / other;
    return out;
}

so& so::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathbb{RR}) \rightarrow \mathfrak{so} \f}
    *
    * In place scalar division.
    */

    this->point /= other;
    return *this;
}

}
