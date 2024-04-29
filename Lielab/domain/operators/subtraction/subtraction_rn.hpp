#ifndef LIELAB_DOMAIN_OPERATORS_subtraction_rn_HPP_
#define LIELAB_DOMAIN_OPERATORS_subtraction_rn_HPP_

#include "../../liealgebras.hpp"
#include "../../liegroups.hpp"

namespace Lielab
{
namespace domain
{

Lielab::domain::cn operator-(const Lielab::domain::rn & lhs, const Lielab::domain::cn & rhs)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{cn}) \rightarrow \mathfrak{glc} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape - 1);
    const Eigen::VectorXcd new_vector = lhs._data(slice).cast<std::complex<double>>() - rhs._data(slice);
    return new_vector;
}

Lielab::domain::gl operator-(const Lielab::domain::rn & lhs, const Lielab::domain::gl & rhs)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{gl}) \rightarrow \mathfrak{gl} \f}
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

Lielab::domain::glc operator-(const Lielab::domain::rn & lhs, const Lielab::domain::glc & rhs)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
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

Lielab::domain::rn operator-(const Lielab::domain::rn & lhs, const Lielab::domain::rn & rhs)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    const size_t new_shape = std::min(lhs.shape, rhs.shape);
    const Eigen::ArithmeticSequence slice = Eigen::seqN(0, new_shape - 1);
    const Eigen::VectorXd new_vector = lhs._data(slice) - rhs._data(slice);
    return new_vector;
}

Lielab::domain::gl operator-(const Lielab::domain::rn & lhs, const Lielab::domain::se & rhs)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{se}) \rightarrow \mathfrak{gl} \f}
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

Lielab::domain::gl operator-(const Lielab::domain::rn & lhs, const Lielab::domain::so & rhs)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{so}) \rightarrow \mathfrak{gl} \f}
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

Lielab::domain::gl operator-(const Lielab::domain::rn & lhs, const Lielab::domain::sp & rhs)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{sp}) \rightarrow \mathfrak{gl} \f}
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

Lielab::domain::glc operator-(const Lielab::domain::rn & lhs, const Lielab::domain::su & rhs)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{su}) \rightarrow \mathfrak{glc} \f}
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
