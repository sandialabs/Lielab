#include "SU.hpp"

#include "SO.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string SU::to_string() const
{
    return "SU(" + std::to_string(this->_shape) + ")";
}

SU::SU() : SU(0)
{
    /*! \f{equation*}{ () \rightarrow SU \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::SU x, y, z;
    * 
    */

}

SU::SU(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SU \f}
    *
    * Constructor instantiating an \f$SU\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::SU x(2), y(4), z(6);
    * 
    * @param[in] shape The shape of the data matrix.
    */

    this->data = Eigen::MatrixXcd::Identity(shape, shape);
    this->_shape = shape;
}

SU SU::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SU \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The SU element. 
    */

    return SU(shape);
}

size_t SU::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    if (this->_shape == 0) return 0;
    return this->_shape * this->_shape - 1;
}

size_t SU::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return this->_shape;
}

size_t SU::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    return static_cast<size_t>(2*std::pow(this->_shape, 2));
}

Eigen::MatrixXcd SU::get_matrix() const
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

    return this->data;
}

SU SU::inverse() const
{
    /*! \f{equation*}{ () \rightarrow SU \f}
    * 
    * Returns the inverse.
    */

    return this->data.inverse();
}

Eigen::VectorXd SU::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */
    
    if (this->_shape == 0) return Eigen::VectorXd::Zero(0);

    const Eigen::MatrixXcd A = this->get_matrix();

    Eigen::VectorXd out = Eigen::VectorXd::Zero(this->get_size());
    size_t kk = 0;
    for (size_t ii = 0; ii < this->_shape; ii++)
    {
        for (size_t jj = 0; jj < this->_shape; jj++)
        {
            out(kk) = std::real(A(ii,jj));
            out(kk+1) = std::imag(A(ii,jj));
            kk = kk + 2;
        }
    }

    return out;
}

void SU::unserialize(const Eigen::VectorXd &vec)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow SU \f}
    * 
    * Sets the SP object from a serialized vector.
    */

    const size_t vdim = vec.size();
    const size_t max_ind = std::min(this->get_size(), vdim);
    
    for (size_t vind = 0; vind < max_ind; vind++)
    {
        const size_t rem = vind % 2;
        const size_t row = static_cast<size_t>(std::floor((vind/2) / this->_shape));
        const size_t col = (vind/2) % this->_shape;
        
        if (rem == 0)
        {
            this->data(row, col).real(vec(vind));
        }
        else
        {
            this->data(row, col).imag(vec(vind));
        }
    }
}

void SU::unserialize(std::initializer_list<double> vec)
{
    this->unserialize(Eigen::VectorXd{std::move(vec)});
}

std::complex<double> SU::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t shape = this->get_shape();
    if (shape == 0) return std::complex<double>(nan, nan);

    if (index1 >= static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);
    if (index2 >= static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);
    if (std::abs(index1) > static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);
    if (std::abs(index2) > static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);

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
    
    return this->data(_index1, _index2);
}

SU SU::operator*(const SU & other) const
{
    /*! \f{equation*}{ (SU, SU) \rightarrow SU \f}
    *
    * Group product.
    */

    assert(this->_shape == other.get_shape());
    return this->data * other.data;
}

SU & SU::operator*=(const SU & other)
{
    /*! \f{equation*}{ (SU, SU) \rightarrow SU \f}
    *
    * In place group product.
    */

    assert(this->_shape == other.get_shape());
    this->data *= other.data;
    return *this;
}

std::ostream & operator<<(std::ostream& os, const SU & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << static_cast<const Eigen::MatrixXcd>(other.data);
    return os;
}

/*
* Additional static initializers. Not a part of the core Lie group, but are convenient.
*/
template SU SU::from_quaternion<double>(const double, const double, const double, const double);

SU SU::from_SO3(const SO & dcm)
{
    /*!
    * Transforms an SO(3) object into an SU(2) object.
    *
    * @param[in] dcm A direction cosine as an SO object of shape 3.
    * @param[out] q A quaternion as an SU object of shape 2.
    */

    if (dcm.get_shape() != 3)
    {
        throw std::domain_error("SU::from_SO3: Expected input shape 3. Got " + std::to_string(dcm.get_shape()) + ".");
    }

    const auto [q0, q1, q2, q3] = dcm.to_quaternion<double>();

    return SU::from_quaternion(q0, q1, q2, q3);

}

std::array<double, 4> SU::to_quaternion() const
{
    /*! \f{equation*}{ (SU) \rightarrow \mathbb{R}^4 \f}
     *
     * Method returning a quaternion as an \f$\mathbb{R}^4\f$ object.
     * 
     * Enables outputs like:
     * 
     *     std::array<double, 4> Quaternion0 = Qobj.to_quaternion();
     * 
     * @param[out] quaternion An array representing the Quaternion.
     */

    if (this->_shape != 2)
    {
        throw std::domain_error("SU::to_quaternion: Expected input shape 2. Got " + std::to_string(this->_shape) + ".");
    }

    const double q0 = (this->data(0,0).real() + this->data(1,1).real())/2.0;
    const double q1 = (this->data(0,0).imag() - this->data(1,1).imag())/2.0;
    const double q2 = (-this->data(0,1).real() + this->data(1,0).real())/2.0;
    const double q3 = (this->data(0,1).imag() + this->data(1,0).imag())/2.0;

    return {q0, q1, q2, q3};
}

}
