#include "RN.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string RN::to_string() const
{
    const size_t shape = this->get_shape();
    if (shape == 0) return "R^nan";
    return "R^" + std::to_string(shape-1);
}

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

RN::RN(const size_t n)
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
    this->data = Eigen::VectorXd::Zero(n);
}

RN RN::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{RN} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
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

// RN::RN(std::initializer_list<double> other)
// {
//     /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow RN \f}
//     *
//     * Constructor instantiating an \f$RN\f$ object from an \f$n \times 1\f$
//     * real vector.
//     *
//     * @param[in] other The object to instantiate from as a real matrix.
//     */

//     this->data = Eigen::VectorXd{std::move(other)};
//     this->_shape = this->data.size() + 1;
// }

Eigen::MatrixXd RN::project(const Eigen::MatrixXd & other)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in RN \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(other.rows(), other.cols());

    Eigen::MatrixXd out = Eigen::MatrixXd::Identity(shape, shape);

    for (size_t ii = 0; ii < shape - 1; ii++)
    {
        out(ii, shape-1) = other(ii, shape-1);
    }

    return out;
}

size_t RN::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    if (this->_shape == 0) return 0; // TODO: Return nan?
    return static_cast<size_t>(this->_shape - 1);
}

size_t RN::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return this->_shape;
}

size_t RN::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    if (this->_shape == 0) return 0;
    return static_cast<size_t>(this->_shape - 1);
}

Eigen::MatrixXd RN::get_matrix() const
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

    for (size_t ii = 0; ii < this->_shape-1; ii++)
    {
        out(ii, this->_shape - 1) = data(ii);
    }

    return out;
}

RN RN::inverse() const
{
    /*! \f{equation*}{ (RN) \rightarrow RN \f}
    * 
    * Returns the inverse.
    */

    return RN::from_vector(-this->data);
}

// Data representation

Eigen::VectorXd RN::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */

    return this->data;
}

void RN::unserialize(const Eigen::VectorXd& vec)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the RN object from a serialized vector.
    */

    const size_t vdim = vec.size();
    const size_t max_ind = std::min(this->get_dimension(), vdim);

    for (size_t vind = 0; vind < max_ind; vind++)
    {
        this->data(vind) = vec(vind);
    }
}

void RN::unserialize(std::initializer_list<double> vec)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->unserialize(Eigen::VectorXd{std::move(vec)});
}

double RN::operator()(const ptrdiff_t index) const
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

    return this->data(_index);
}

double RN::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t shape = this->get_shape();
    if (shape == 0) return nan;

    if (index1 >= static_cast<ptrdiff_t>(shape)) return nan;
    if (index2 >= static_cast<ptrdiff_t>(shape)) return nan;
    if (std::abs(index1) > static_cast<ptrdiff_t>(shape)) return nan;
    if (std::abs(index2) > static_cast<ptrdiff_t>(shape)) return nan;

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

    if (_index1 == _index2) return 1.0;

    if (_index1 == shape - 1) return 0.0;
    if (_index2 != shape - 1) return 0.0;

    return this->data(_index1);
}

RN RN::operator*(const RN & other) const
{
    /*! \f{equation*}{ (RN, RN) \rightarrow RN \f}
    *
    * Group product.
    */

    assert(this->_shape == other.get_shape());
    return RN::from_vector(this->data + other.data);
}

RN & RN::operator*=(const RN & other)
{
    /*! \f{equation*}{ (RN, RN) \rightarrow RN \f}
    *
    * In place group product.
    */

    assert(this->_shape == other.get_shape());
    this->data += other.data;
    return *this;
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

    const size_t n = other.size();
    RN out(n);
    out.data = other;
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

std::ostream & operator<<(std::ostream& os, const RN & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << static_cast<const Eigen::VectorXd>(other.data);
    return os;
}

}
