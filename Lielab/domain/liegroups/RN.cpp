#include "RN.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>


namespace Lielab::domain
{

RN::RN() : RN(0)
{
    /*! \f{equation*}{ () \rightarrow RN \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::RN x, y, z;
    * 
    */

}

// RN::~RN();

RN::RN(const RN::matrix_t& matrix)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& RN \\ (\mathbb{R}^{n \times 1}) &\rightarrow& RN \f}
    *
    * Constructor instantiating an \f$RN\f$ object from either an
    * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->_shape = static_cast<int>(matrix.rows());

    if (this->_shape == 0)
    {
        this->point.noalias() = Eigen::VectorXd::Zero(0);
        return;
    }

    this->point = Eigen::VectorXd::Zero(this->_shape - 1);
    for (int ii = 0; ii < this->_shape - 1; ii++)
    {
        this->point(ii) = matrix(ii, this->_shape - 1);
    }
}

RN RN::identity(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{RN} \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The RN element. 
    */

    if (shape == 0)
    {
        RN out;
        out._shape = 0;
        return out;
    }

    return RN(shape - 1);
}

RN RN::project(const RN::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in RN \f}
    */

    const size_t shape = std::min(matrix.rows(), matrix.cols());

    Eigen::MatrixXd out = Eigen::MatrixXd::Identity(shape, shape);

    for (size_t ii = 0; ii < shape - 1; ii++)
    {
        out(ii, shape-1) = matrix(ii, shape-1);
    }

    return RN(out);
}

RN::RN(const int n)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow RN \f}
    *
    * Constructor instantiating an \f$RN\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::RN x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */
    
    this->_shape = n + 1;
    this->point.noalias() = Eigen::VectorXd::Zero(n);
}

RN RN::from_vector(const Eigen::VectorXd& other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{rn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{rn}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    const int n = static_cast<int>(other.size());
    RN out(n);
    out.point = other;
    return out;
}

RN RN::from_vector(const std::initializer_list<double> other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{rn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{rn}\f$ object from a
    * \f$n \times 1\f$ vector.
    *
    * @param[in] other The object to instantiate from as a vector.
    */

    return RN::from_vector(Eigen::VectorXd{std::move(other)});
}

std::string RN::to_string() const
{
    const size_t shape = this->get_shape();
    if (shape == 0) return "R^nan";
    return "R^" + std::to_string(shape-1);
}

int RN::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    return static_cast<int>(this->point.size());
}

int RN::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    return static_cast<int>(this->point.size());
}

bool RN::is_abelian() const
{
    return true;
}

int RN::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return this->_shape;
}

RN::point_t RN::get_point() const
{
    return this->point;
}

Eigen::VectorXd RN::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */

    return this->point;
}

void RN::unserialize(const Eigen::VectorXd& serialized)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the RN object from a serialized vector.
    */

    const int vdim = static_cast<int>(serialized.size());
    const int max_ind = std::min(this->get_dimension(), vdim);

    for (int vind = 0; vind < max_ind; vind++)
    {
        this->point(vind) = serialized(vind);
    }
}

void RN::unserialize(std::initializer_list<double> serialized)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

RN::matrix_t RN::get_matrix() const
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

    Eigen::MatrixXd out = Eigen::MatrixXd::Identity(this->_shape, this->_shape);

    if (this->_shape == 0) return out;

    for (int ii = 0; ii < this->_shape-1; ii++)
    {
        out(ii, this->_shape - 1) = this->point(ii);
    }

    return out;
}

RN::field_t RN::operator()(const int index1, const int index2) const
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

    if (_index1 == _index2) return 1.0;
    if (_index1 == shape - 1) return 0.0;
    if (_index2 != shape - 1) return 0.0;

    return this->point(_index1);
}

const RN::field_t& RN::operator[](const int index) const
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    */

    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < dim), "Index " + std::to_string(index) + " is out of bounds for RN of length " + std::to_string(dim));

    return this->point(_index);
}

RN::field_t& RN::operator[](const int index)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    */

    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < dim), "Index " + std::to_string(index) + " is out of bounds for RN of length " + std::to_string(dim));

    return this->point(_index);
}

RN RN::operator*(const RN& other) const
{
    /*! \f{equation*}{ (RN, RN) \rightarrow RN \f}
    *
    * Group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return RN::from_vector(this->point + other.point);
}

RN& RN::operator*=(const RN& other)
{
    /*! \f{equation*}{ (RN, RN) \rightarrow RN \f}
    *
    * In place group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point += other.point;
    return *this;
}

RN RN::inverse() const
{
    /*! \f{equation*}{ (RN) \rightarrow RN \f}
    * 
    * Returns the inverse.
    */

    return RN::from_vector(-this->point);
}

}
