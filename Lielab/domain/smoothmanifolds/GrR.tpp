#ifndef LIELAB_DOMAIN_GRR_TPP
#define LIELAB_DOMAIN_GRR_TPP

#include "GrR.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>

namespace Lielab::domain
{

template<typename OtherDerived>
GrR::GrR(const Eigen::MatrixBase<OtherDerived>& other)
{
    /*! 
    */

    if (other.cols() > other.rows())
    {
        throw Lielab::utils::Error("n must be greater than or equal to k.");
    }

    this->data = Eigen::MatrixXd(other);
}

template<typename OtherDerived>
GrR& GrR::operator=(const Eigen::MatrixBase<OtherDerived>& other)
{
    /*! 
    */

    if (other.cols() != other.rows())
    {
        throw Lielab::utils::Error("n must be greater than or equal to k.");
    }
    
    this->data = Eigen::MatrixXd(other);
    return *this;
}

}

#endif
