#ifndef LIELAB_DOMAIN_OPERATORS_alg_grp_prod_sp_HPP_
#define LIELAB_DOMAIN_OPERATORS_alg_grp_prod_sp_HPP_

#include "../../liealgebras.hpp"
#include "../../liegroups.hpp"

namespace Lielab
{
namespace domain
{

Lielab::domain::glc operator*(const Lielab::domain::sp & lhs, const Lielab::domain::CN & rhs)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{CN}) \rightarrow \mathfrak{glc} \f}
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

Lielab::domain::gl operator*(const Lielab::domain::sp & lhs, const Lielab::domain::GL & rhs)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{GL}) \rightarrow \mathfrak{gl} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator*(const Lielab::domain::sp & lhs, const Lielab::domain::GLC & rhs)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{GLC}) \rightarrow \mathfrak{glc} \f}
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

Lielab::domain::gl operator*(const Lielab::domain::sp & lhs, const Lielab::domain::RN & rhs)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{RN}) \rightarrow \mathfrak{gl} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::gl operator*(const Lielab::domain::sp & lhs, const Lielab::domain::SE & rhs)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{SE}) \rightarrow \mathfrak{gl} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::gl operator*(const Lielab::domain::sp & lhs, const Lielab::domain::SO & rhs)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{SO}) \rightarrow \mathfrak{gl} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::sp operator*(const Lielab::domain::sp & lhs, const Lielab::domain::SP & rhs)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{SP}) \rightarrow \mathfrak{sp} \f}
    *
    * Vector product.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) * rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator*(const Lielab::domain::sp & lhs, const Lielab::domain::SU & rhs)
{
    /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{SU}) \rightarrow \mathfrak{glc} \f}
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

}
}

#endif
