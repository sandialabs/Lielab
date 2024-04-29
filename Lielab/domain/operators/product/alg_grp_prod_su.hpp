#ifndef LIELAB_DOMAIN_OPERATORS_alg_grp_prod_su_HPP_
#define LIELAB_DOMAIN_OPERATORS_alg_grp_prod_su_HPP_

#include "../../liealgebras.hpp"
#include "../../liegroups.hpp"

namespace Lielab
{
namespace domain
{

Lielab::domain::glc operator*(const Lielab::domain::su & lhs, const Lielab::domain::CN & rhs)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{CN}) \rightarrow \mathfrak{glc} \f}
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

Lielab::domain::glc operator*(const Lielab::domain::su & lhs, const Lielab::domain::GL & rhs)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{GL}) \rightarrow \mathfrak{glc} \f}
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

Lielab::domain::glc operator*(const Lielab::domain::su & lhs, const Lielab::domain::GLC & rhs)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{GLC}) \rightarrow \mathfrak{glc} \f}
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

Lielab::domain::glc operator*(const Lielab::domain::su & lhs, const Lielab::domain::RN & rhs)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{RN}) \rightarrow \mathfrak{glc} \f}
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

Lielab::domain::glc operator*(const Lielab::domain::su & lhs, const Lielab::domain::SE & rhs)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{SE}) \rightarrow \mathfrak{glc} \f}
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

Lielab::domain::glc operator*(const Lielab::domain::su & lhs, const Lielab::domain::SO & rhs)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{SO}) \rightarrow \mathfrak{glc} \f}
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

Lielab::domain::glc operator*(const Lielab::domain::su & lhs, const Lielab::domain::SP & rhs)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{SP}) \rightarrow \mathfrak{glc} \f}
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

Lielab::domain::su operator*(const Lielab::domain::su & lhs, const Lielab::domain::SU & rhs)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{SU}) \rightarrow \mathfrak{su} \f}
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
