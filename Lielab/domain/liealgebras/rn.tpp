#ifndef LIELAB_DOMAIN_LIEALGEBRAS_rn_TPP
#define LIELAB_DOMAIN_LIEALGEBRAS_rn_TPP

#include "rn.hpp"

#include "LieAlgebra.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
rn::rn(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{rn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{rn}\f$ object from an
    * \f$n \times n\f$ real matrix.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Input data to rn malformed.");
    }

    this->_shape = other.rows();

    if (this->_shape == 0)
    {
        this->data = Eigen::VectorXd::Zero(0);
        return;
    }

    this->data = Eigen::VectorXd::Zero(this->_shape - 1);
    for (size_t ii = 0; ii < this->_shape - 1; ii++)
    {
        this->data(ii) = other(ii, this->_shape - 1);
    }
}

}

#endif
