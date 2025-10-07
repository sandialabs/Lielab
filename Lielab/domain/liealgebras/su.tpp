#ifndef LIELAB_DOMAIN_LIEALGEBRAS_su_TPP
#define LIELAB_DOMAIN_LIEALGEBRAS_su_TPP

#include "su.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
su::su(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{equation}{ (\mathbb{C}^{n \times n}) \rightarrow \mathfrak{su} \f}
    *
    * Constructor instantiating an \f$\mathfrak{su}\f$ object from an
    * \f$n \times n\f$ real matrix.
    *
    * @param[in] other The object to instantiate from.
    */
    
    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Input data to su malformed.");
    }

    this->data = Eigen::MatrixXcd(other);
    this->_shape = this->data.rows();
}

}

#endif
