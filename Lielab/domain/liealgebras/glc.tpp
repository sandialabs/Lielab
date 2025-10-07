#ifndef LIELAB_DOMAIN_LIEALGEBRAS_glc_TPP
#define LIELAB_DOMAIN_LIEALGEBRAS_glc_TPP

#include "glc.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
glc::glc(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{eqnarray*}{(\mathbb{C}^{n \times n}) &\rightarrow& \mathfrak{glc} \f}
    *
    * Constructor instantiating an \f$\mathfrak{glc}\f$ object from an
    * \f$n \times n\f$ imaginary matrix.
    *
    * @param[in] other The object to instantiate from as an imaginary matrix.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Input data to glc malformed.");
    }

    this->_shape = other.rows();
    this->data = other;
}

}

#endif
