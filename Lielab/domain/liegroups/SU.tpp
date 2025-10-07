#ifndef LIELAB_DOMAIN_SU_TPP
#define LIELAB_DOMAIN_SU_TPP

#include "SU.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
SU::SU(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{equation}{(\mathbb{C}^{n \times n}) \rightarrow SU \f}
    *
    * Constructor instantiating an \f$SU\f$ object from an
    * \f$n \times n\f$ imaginary matrix.
    *
    * @param[in] other The object to instantiate from as an imaginary matrix.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Size of the matrix must be square.");
    }

    this->data = Eigen::MatrixXcd(other);
    this->_shape = other.rows();
}

template<typename OtherDerived>
SU & SU::operator=(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{equation}{ SU := \mathbb{C}^{n \times n} \f}
    * 
    * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `SU`.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Size of the matrix must be square.");
    }
    
    this->data = Eigen::MatrixXcd(other);
    this->_shape = other.rows();
    return *this;
}

/*
* Additional static initializers. Not a part of the core Lie group, but are convenient.
*/
template <typename T>
SU SU::from_quaternion(const T e0, const T e1, const T e2, const T e3)
{
    /*! \f{equation*}{ (\mathbb{R}^4) \rightarrow SU \f}
     *
     * Constructor instantiating a Quaternion as an \f$SU\f$ object.
     * 
     * Enables instatiation like:
     * 
     *     Lielab::domain::SU Quaternion0 = Lielab::domain::SU::from_quaternion(1.0, 0.0, 0.0, 0.0);
     * 
     * @param[out] quaternion An SU object representing the Quaternion.
     */

    constexpr std::complex<T> j(0.0, 1.0);
    SU qout = SU(2);

    qout.data(0,0) = e0 + e1*j;
    qout.data(1,1) = e0 - e1*j;
    qout.data(0,1) = -e2 + e3*j;
    qout.data(1,0) = e2 + e3*j;

    return qout;
}

}

#endif
