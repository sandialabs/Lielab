#ifndef LIELAB_DOMAIN_GL_TPP
#define LIELAB_DOMAIN_GL_TPP

#include "GLR.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
GLR::GLR(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& GLR \\ (\mathbb{R}^{n \times 1}) &\rightarrow& GLR \f}
    *
    * Constructor instantiating an \f$GLR\f$ object from either an
    * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
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
        throw Lielab::utils::Error("Input data to GLR malformed.");
    }
}

template<typename OtherDerived>
GLR & GLR::operator=(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{eqnarray*}{ GLR &:= \mathbb{R}^{n \times n} \\ &:= \mathbb{R}^{n \times 1} \f}
    * 
    * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `GLR`.
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
