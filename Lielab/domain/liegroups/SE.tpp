#ifndef LIELAB_DOMAIN_SE_TPP
#define LIELAB_DOMAIN_SE_TPP

#include "SE.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

template<typename OtherDerived>
SE::SE(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{equation}{ (\mathbb{R}^{n \times n}) \rightarrow SE \f}
    *
    * Constructor instantiating an \f$SE\f$ object from an
    * \f$n \times n\f$ real matrix.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Size of the matrix must be square.");
    }

    this->data = Eigen::MatrixXd(other);
    this->_shape = this->data.rows();
}

template<typename OtherDerived>
SE & SE::operator=(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{equation}{ SE := \mathbb{R}^{n \times n} \f}
    * 
    * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `SE`.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Size of the matrix must be square.");
    }
    
    this->data = data_t(other);
    this->_shape = this->data.rows();
    return *this;
}

}

#endif
