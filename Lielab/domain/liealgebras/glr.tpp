#ifndef LIELAB_DOMAIN_LIEALGEBRAS_glr_TPP
#define LIELAB_DOMAIN_LIEALGEBRAS_glr_TPP

#include "glr.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
glr::glr(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{glr} \\ (\mathbb{R}^{n \times 1}) &\rightarrow& \mathfrak{glr} \f}
    *
    * Constructor instantiating an \f$\mathfrak{glr}\f$ object from either an
    * \f$n \times n\f$ real matrix.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Input data to glr malformed.");
    }

    this->_shape = other.rows();
    this->data = other;
}

}

#endif
