#ifndef LIELAB_DOMAIN_CN_TPP
#define LIELAB_DOMAIN_CN_TPP

#include "CN.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
CN::CN(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{eqnarray*}{(\mathbb{C}^{n \times n}) &\rightarrow& CN \\ (\mathbb{C}^{n \times 1}) &\rightarrow& CN \f}
    *
    * Constructor instantiating an \f$CN\f$ object from either an
    * \f$n \times n\f$ imaginary matrix or \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Input data to CN malformed.");
    }

    this->_shape = other.rows();

    if (this->_shape == 0)
    {
        this->data = Eigen::VectorXcd::Zero(0);
        return;
    }

    this->data = Eigen::VectorXcd::Zero(this->_shape - 1);
    for (size_t ii = 0; ii < this->_shape - 1; ii++)
    {
        this->data(ii) = other(ii, this->_shape - 1);
    }
}

}

#endif
