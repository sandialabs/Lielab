#ifndef LIELAB_DOMAIN_LIEALGEBRAS_sp_TPP
#define LIELAB_DOMAIN_LIEALGEBRAS_sp_TPP

#include "sp.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
sp::sp(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{equation}{ (\mathbb{R}^{n \times n}) \rightarrow \mathfrak{sp} \f}
    *
    * Constructor instantiating an \f$\mathfrak{sp}\f$ object from an
    * \f$n \times n\f$ real matrix.
    *
    * @param[in] other The object to instantiate from.
    */
    
    if (other.rows() % 2 != 0)
    {
        throw Lielab::utils::Error("Shape of sp must be even dimensional.");
    }

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Shape of sp must be square.");
    }

    this->data = Eigen::MatrixXd(other);
    this->_shape = this->data.rows();
}

}

#endif
