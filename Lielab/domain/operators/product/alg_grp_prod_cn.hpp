#ifndef LIELAB_DOMAIN_OPERATORS_alg_grp_prod_cn_HPP_
#define LIELAB_DOMAIN_OPERATORS_alg_grp_prod_cn_HPP_

#include "../../liealgebras.hpp"
#include "../../liegroups.hpp"

namespace Lielab
{
namespace domain
{

Lielab::domain::glc operator*(const Lielab::domain::cn & lhs, const Lielab::domain::CN & rhs)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{CN}) \rightarrow \mathfrak{glc} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator*(const Lielab::domain::cn & lhs, const Lielab::domain::GL & rhs)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{GL}) \rightarrow \mathfrak{glc} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator*(const Lielab::domain::cn & lhs, const Lielab::domain::GLC & rhs)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{GLC}) \rightarrow \mathfrak{glc} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator*(const Lielab::domain::cn & lhs, const Lielab::domain::RN & rhs)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{RN}) \rightarrow \mathfrak{glc} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator*(const Lielab::domain::cn & lhs, const Lielab::domain::SE & rhs)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{SE}) \rightarrow \mathfrak{glc} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator*(const Lielab::domain::cn & lhs, const Lielab::domain::SO & rhs)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{SO}) \rightarrow \mathfrak{glc} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator*(const Lielab::domain::cn & lhs, const Lielab::domain::SP & rhs)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{SP}) \rightarrow \mathfrak{glc} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator*(const Lielab::domain::cn & lhs, const Lielab::domain::SU & rhs)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{SU}) \rightarrow \mathfrak{glc} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

}
}

#endif
