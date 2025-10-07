#ifndef LIELAB_DOMAIN_LIEALGEBRAS_se_TPP
#define LIELAB_DOMAIN_LIEALGEBRAS_se_TPP

#include "se.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
se::se(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{se} \f}
    *
    * Constructor instantiating an \f$\mathfrak{se}\f$ object from either an
    * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Input data to se malformed.");
    }
    
    this->data = Eigen::MatrixXd(other);
    this->_shape = this->data.rows();
}

}

#endif
