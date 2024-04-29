#ifndef LIELAB_DOMAIN_OPERATORS_addition_glc_HPP_
#define LIELAB_DOMAIN_OPERATORS_addition_glc_HPP_

#include "../../liealgebras.hpp"
#include "../../liegroups.hpp"

namespace Lielab
{
namespace domain
{

Lielab::domain::glc operator+(const Lielab::domain::glc & lhs, const Lielab::domain::cn & rhs)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{cn}) \rightarrow \mathfrak{glc} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator+(const Lielab::domain::glc & lhs, const Lielab::domain::gl & rhs)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{gl}) \rightarrow \mathfrak{glc} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator+(const Lielab::domain::glc & lhs, const Lielab::domain::glc & rhs)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator+(const Lielab::domain::glc & lhs, const Lielab::domain::rn & rhs)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{rn}) \rightarrow \mathfrak{glc} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator+(const Lielab::domain::glc & lhs, const Lielab::domain::se & rhs)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{se}) \rightarrow \mathfrak{glc} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator+(const Lielab::domain::glc & lhs, const Lielab::domain::so & rhs)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{so}) \rightarrow \mathfrak{glc} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator+(const Lielab::domain::glc & lhs, const Lielab::domain::sp & rhs)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{sp}) \rightarrow \mathfrak{glc} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix().cast<std::complex<double>>();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return new_matrix;
}

Lielab::domain::glc operator+(const Lielab::domain::glc & lhs, const Lielab::domain::su & rhs)
{
    /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{so}) \rightarrow \mathfrak{glc} \f}
    *
    * Addition of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape);
    const Eigen::MatrixXcd lhs_matrix = lhs.get_matrix();
    const Eigen::MatrixXcd rhs_matrix = rhs.get_matrix();
    const Eigen::MatrixXcd new_matrix = lhs_matrix(slice, slice) + rhs_matrix(slice, slice);

    return new_matrix;
}

}
}

#endif
