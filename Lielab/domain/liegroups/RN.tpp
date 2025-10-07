#ifndef LIELAB_DOMAIN_RN_TPP
#define LIELAB_DOMAIN_RN_TPP

#include "RN.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
RN::RN(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& RN \\ (\mathbb{R}^{n \times 1}) &\rightarrow& RN \f}
    *
    * Constructor instantiating an \f$RN\f$ object from either an
    * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Input data to RN malformed.");
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
