#ifndef LIELAB_DOMAIN_SP_TPP
#define LIELAB_DOMAIN_SP_TPP

#include "SP.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
SP::SP(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times n}) \rightarrow SP \f}
    *
    * Constructor instantiating an \f$SP\f$ object from an
    * \f$n \times n\f$ real matrix.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Size of the matrix must be square.");
    }

    if (other.rows() % 2 != 0)
    {
        throw Lielab::utils::Error("Input matrix must be even dimensional.");
    }

    this->data = Eigen::MatrixXd(other);
    this->_shape = other.rows();
}

template<typename OtherDerived>
SP & SP::operator=(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{equation}{ SP := \mathbb{R}^{n \times n} \f}
    * 
    * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `SP`.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Size of the matrix must be square.");
    }

    if (other.rows() % 2 != 0)
    {
        throw Lielab::utils::Error("Input matrix must be even dimensional.");
    }
    
    this->data = Eigen::MatrixXd(other);
    this->_shape = this->data.rows();
    return *this;
}

}

#endif
