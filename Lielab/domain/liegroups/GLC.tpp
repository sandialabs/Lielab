#ifndef LIELAB_DOMAIN_GLC_TPP
#define LIELAB_DOMAIN_GLC_TPP

#include "GLC.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
GLC::GLC(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& GLC \\ (\mathbb{R}^{n \times 1}) &\rightarrow& GLC \f}
    *
    * Constructor instantiating an \f$GLC\f$ object from either an
    * \f$n \times n\f$ imaginary matrix or \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    if (other.rows() == other.cols())
    {
        this->data = other;
        this->_shape = other.rows();
    }
    else
    {
        throw Lielab::utils::Error("Input data to GLC malformed.");
    }
}

template<typename OtherDerived>
GLC & GLC::operator=(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{eqnarray*}{ GLC &:= \mathbb{C}^{n \times n} \\ &:= \mathbb{C}^{n \times 1} \f}
    * 
    * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `GLC`.
    */

    if (other.rows() == other.cols())
    {   
        this->_shape = other.rows();
        this->data = other;
    }
    return *this;
}

}

#endif
