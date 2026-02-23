#include "rn.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>


namespace Lielab::domain
{

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

// rn::~rn();

rn::rn(const rn::matrix_t& matrix)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{rn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{rn}\f$ object from an
    * \f$n \times n\f$ real matrix.
    *
    * @param[in] matrix The object to instantiate from as a real matrix.
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->_shape = static_cast<int>(matrix.rows());

    if (this->_shape == 0)
    {
        this->point.noalias() = Eigen::VectorXd::Zero(0);
        return;
    }

    this->point.noalias() = Eigen::VectorXd(this->_shape - 1);
    for (int ii = 0; ii < this->_shape - 1; ii++)
    {
        this->point(ii) = matrix(ii, this->_shape - 1);
    }
}

rn rn::basis(const int index, const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{rn} \f}
    *
    * Returns the i'th basis element of the rn algebra.
    * 
    * @param[in] index The basis vector.
    * @param[in] shape The shape of the algebra.
    * @param[out] out The rn element. 
    */

    rn out(shape);
    if (index < 0) return out;

    const int dim = out.get_dimension();
    if (index >= dim) return out;

    out.point(index) = 1.0;

    return out;
}

rn rn::zero(const int shape)
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

rn rn::from_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{rn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{rn}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    const int shape = static_cast<int>(vector.size()) + 1;
    rn out = rn::zero(shape);
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

rn rn::project(const rn::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{rn} \f}
    *
    */

    const size_t shape = std::min(matrix.rows(), matrix.cols());
    return rn(matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape)));
}

rn::rn(const int n)
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
    this->point.noalias() = Eigen::VectorXd::Zero(n);
}

std::string rn::to_string() const
{
    if (this->_shape == 0) return "r^nan";
    return "r^" + std::to_string(this->_shape-1);
}

int rn::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    if (this->_shape == 0) return 0; // Return nan?

    return this->_shape - 1;
}

int rn::get_size() const
{
    return this->get_dimension();
}

bool rn::is_abelian() const
{
    return true;
}

int rn::get_shape() const
{
    return this->_shape;
}

rn::point_t rn::get_point() const
{
    /*!
    */

    return this->point;
}

Eigen::VectorXd rn::serialize() const
{
    return this->get_vector();
}

void rn::unserialize(const Eigen::VectorXd& serialized)
{
    this->set_vector(serialized);
}

void rn::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
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

    for (int ii = 0; ii < this->_shape-1; ii++)
    {
        out(ii, this->_shape - 1) = this->point(ii);
    }

    return out;
}

Eigen::VectorXd rn::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */

    return this->point;
}

void rn::set_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ \mathfrak{rn} := \mathbb{R}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    */

    const int vdim = static_cast<int>(vector.size());
    const int max_ind = std::min(this->get_dimension(), vdim);

    for (int vind = 0; vind < max_ind; vind++)
    {
        this->point(vind) = vector(vind);
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

double rn::operator()(const int index) const
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

    return this->point(_index);
}

rn::field_t rn::operator()(const int index1, const int index2) const
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

    if (_index1 == shape - 1) return 0.0;
    if (_index2 != shape - 1) return 0.0;

    return this->point(_index1);
}

const rn::field_t& rn::operator[](const int index) const
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    */

    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < dim), "Index " + std::to_string(index) + " is out of bounds for rn of length " + std::to_string(dim));

    return this->point(_index);
}

rn::field_t& rn::operator[](const int index)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    */

    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < dim), "Index " + std::to_string(index) + " is out of bounds for rn of length " + std::to_string(dim));

    return this->point(_index);
}

rn rn::operator+(const rn& other) const
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * Addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return rn::from_vector(this->point + other.point);
}

rn& rn::operator+=(const rn& other)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point += other.point;
    return *this;
}

rn rn::operator-(const rn& other) const
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return rn::from_vector(this->point - other.point);
}

rn& rn::operator-=(const rn& other)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point -= other.point;
    return *this;
}

rn rn::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * Unary negative of the vector.
    */

    return rn::from_vector(-this->point);
}

rn rn::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
    *
    * Scalar product.
    */

    return rn::from_vector(this->point * other);
}

rn operator*(const double other, const rn& rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * Scalar product.
    */

    return rn::from_vector(rhs.point * other);
}

rn& rn::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
    *
    * In place scalar product.
    */

    this->point *= other;
    return *this;
}

rn rn::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
    *
    * Scalar division.
    */

    return rn::from_vector(this->point / other);
}

rn& rn::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
    *
    * In place scalar division.
    */

    this->point /= other;
    return *this;
}

}
