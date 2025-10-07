#ifndef LIELAB_DOMAIN_LIEALGEBRAS_cn_TPP
#define LIELAB_DOMAIN_LIEALGEBRAS_cn_TPP

#include "LieAlgebra.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <complex>

namespace Lielab::domain
{

template<typename OtherDerived>
cn::cn(const Eigen::MatrixBase<OtherDerived>& other)
{
    /*! \f{eqnarray*}{(\mathbb{C}^{n \times n}) &\rightarrow& \mathfrak{cn} \\ (\mathbb{C}^{n \times 1}) &\rightarrow& \mathfrak{cn} \f}
    *
    * Constructor instantiating an \f$\mathfrak{cn}\f$ object from either an
    * \f$n \times n\f$ imaginary matrix or \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary matrix.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Input data to cn malformed.");
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
