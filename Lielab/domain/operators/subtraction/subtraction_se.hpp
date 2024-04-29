#ifndef LIELAB_DOMAIN_OPERATORS_subtraction_se_HPP_
#define LIELAB_DOMAIN_OPERATORS_subtraction_se_HPP_

#include "../../liealgebras.hpp"
#include "../../liegroups.hpp"

namespace Lielab
{
namespace domain
{

Lielab::domain::glc operator-(const Lielab::domain::se & lhs, const Lielab::domain::cn & rhs)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{cn}) \rightarrow \mathfrak{glc} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::gl operator-(const Lielab::domain::se & lhs, const Lielab::domain::gl & rhs)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{gl}) \rightarrow \mathfrak{gl} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator-(const Lielab::domain::se & lhs, const Lielab::domain::glc & rhs)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::gl operator-(const Lielab::domain::se & lhs, const Lielab::domain::rn & rhs)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{rn}) \rightarrow \mathfrak{gl} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::se operator-(const Lielab::domain::se & lhs, const Lielab::domain::se & rhs)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::gl operator-(const Lielab::domain::se & lhs, const Lielab::domain::so & rhs)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{so}) \rightarrow \mathfrak{gl} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::gl operator-(const Lielab::domain::se & lhs, const Lielab::domain::sp & rhs)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{sp}) \rightarrow \mathfrak{gl} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator-(const Lielab::domain::se & lhs, const Lielab::domain::su & rhs)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{su}) \rightarrow \mathfrak{glc} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) - rhs_matrix(slice, slice);

    return new_matrix;
}

}
}

#endif
