#include "SO.hpp"

#include "SU.hpp"
#include "Lielab/domain/liealgebras/so.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <cassert>

namespace Lielab::domain
{

SO::SO() : SO(0)
{
    /*! \f{equation*}{ () \rightarrow SO \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::SO x, y, z;
    * 
    */

}

// SO::~SO();

SO::SO(const SO::matrix_t& matrix)
{
    /*! \f{equation}{ (\mathbb{R}^{n \times n}) \rightarrow SO \f}
    *
    * Constructor instantiating an \f$SO\f$ object from an
    * \f$n \times n\f$ real matrix.
    *
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->point.noalias() = matrix;
}

SO SO::identity(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SO \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The SO element. 
    */

    return SO(shape);
}

SO SO::project(const SO::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in SO \f}
    *
    */

    const int shape = static_cast<int>(std::min(matrix.rows(), matrix.cols()));

    if (shape == 0) return SO::identity(0);

    const Eigen::MatrixXd matrix_sq = matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape));

    Eigen::JacobiSVD<SO::matrix_t, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(matrix_sq);
    SO::matrix_t Q = svd.matrixU() * svd.matrixV().transpose();

    if (Q.determinant() < 0.0)
    {
        int min_col = 0;
        double min_norm = Q.col(0).norm();

        for (int ii = 1; ii < shape; ii++)
        {   
            const double _norm = Q.col(ii).norm();
            if (_norm < min_norm)
            {
                min_col = ii;
                min_norm = _norm;
            }
        }

        Q.col(min_col) = -Q.col(min_col);
    }

    return SO(Q);
}

SO::SO(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SO \f}
    *
    * Constructor instantiating an \f$SO\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::SO x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */

    this->point.noalias() = point_t::Identity(shape, shape);
}

SO SO::from_eulerangles_body123(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_2 \cos \theta_3 & -\cos \theta_2 \sin \theta_3 & \sin \theta_2 \\
    * \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_2 \\
    * \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_2
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-123 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    *
    * TODO: Citation needed.
    *
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta2*ctheta3;
    dcm.point(0,1) = -ctheta2*stheta3;
    dcm.point(0,2) = stheta2;
    dcm.point(1,0) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(1,1) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
    dcm.point(1,2) = -stheta1*ctheta2;
    dcm.point(2,0) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.point(2,1) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(2,2) = ctheta1*ctheta2;

    return dcm;
}


SO SO::from_eulerangles_body231(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_2 & \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 \\
    * \sin \theta_2 & \cos \theta_2 \cos \theta_3 & - \cos \theta_2 \sin \theta_3 \\
    * - \sin \theta_1 \cos \theta_2 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-231 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta2;
    dcm.point(0,1) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.point(0,2) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(1,0) = stheta2;
    dcm.point(1,1) = ctheta2*ctheta3;
    dcm.point(1,2) = -ctheta2*stheta3;
    dcm.point(2,0) = -stheta1*ctheta2;
    dcm.point(2,1) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(2,2) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;

    return dcm;
}


SO SO::from_eulerangles_body312(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_2 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3 \\
    * \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_2 & \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 \\
    * - \cos \theta_2 \sin \theta_3 & \sin \theta_2 & \cos \theta_2 \cos \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-312 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
    dcm.point(0,1) = -stheta1*ctheta2;
    dcm.point(0,2) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(1,0) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(1,1) = ctheta1*ctheta2;
    dcm.point(1,2) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.point(2,0) = -ctheta2*stheta3;
    dcm.point(2,1) = stheta2;
    dcm.point(2,2) = ctheta2*ctheta3;

    return dcm;
}


SO SO::from_eulerangles_body132(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_2 \cos \theta_3 & - \sin \theta_2 & \cos \theta_2 \sin \theta_3 \\
    * \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_2 & - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 \\
    * - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \cos \theta_2 & \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-132 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta2*ctheta3;
    dcm.point(0,1) = -stheta2;
    dcm.point(0,2) = ctheta2*stheta3;
    dcm.point(1,0) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.point(1,1) = ctheta1*ctheta2;
    dcm.point(1,2) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(2,0) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(2,1) = stheta1*ctheta2;
    dcm.point(2,2) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;

    return dcm;
}


SO SO::from_eulerangles_body213(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \cos \theta_2 \\
    * \cos \theta_2 \sin \theta_3 & \cos \theta_2 \cos \theta_3 & - \sin \theta_2 \\
    * - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_2
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-213 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
    dcm.point(0,1) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(0,2) = stheta1*ctheta2;
    dcm.point(1,0) = ctheta2*stheta3;
    dcm.point(1,1) = ctheta2*ctheta3;
    dcm.point(1,2) = -stheta2;
    dcm.point(2,0) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(2,1) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.point(2,2) = ctheta1*ctheta2;

    return dcm;
}


SO SO::from_eulerangles_body321(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_2 & - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 \\
    * \sin \theta_1 \cos \theta_2 & \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 \\
    * - \sin \theta_2 & \cos \theta_2 \sin \theta_3 & \cos \theta_2 \cos \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-321 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta2;
    dcm.point(0,1) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(0,2) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.point(1,0) = stheta1*ctheta2;
    dcm.point(1,1) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
    dcm.point(1,2) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(2,0) = -stheta2;
    dcm.point(2,1) = ctheta2*stheta3;
    dcm.point(2,2) = ctheta2*ctheta3;

    return dcm;
}


SO SO::from_eulerangles_body121(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_2 & \sin \theta_2 \sin \theta_3 & \sin \theta_2 \cos \theta_3 \\
    * \sin \theta_1 \sin \theta_2 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 \\
    * - \cos \theta_1 \sin \theta_2 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-121 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta2;
    dcm.point(0,1) = stheta2*stheta3;
    dcm.point(0,2) = stheta2*ctheta3;
    dcm.point(1,0) = stheta1*stheta2;
    dcm.point(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.point(1,2) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.point(2,0) = -ctheta1*stheta2;
    dcm.point(2,1) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.point(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

    return dcm;
}


SO SO::from_eulerangles_body131(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_2 & - \sin \theta_2 \cos \theta_3 & \sin \theta_2 \sin \theta_3 \\
    * \cos \theta_1 \sin \theta_2 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 \\
    * \sin \theta_1 \sin \theta_2 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-131 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta2;
    dcm.point(0,1) = -stheta2*ctheta3;
    dcm.point(0,2) = stheta2*stheta3;
    dcm.point(1,0) = ctheta1*stheta2;
    dcm.point(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.point(1,2) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.point(2,0) = stheta1*stheta2;
    dcm.point(2,1) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.point(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

    return dcm;
}


SO SO::from_eulerangles_body212(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_1 \sin \theta_2 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 \\
    * \sin \theta_2 \sin \theta_3 & \cos \theta_2 & - \sin \theta_2 \cos \theta_3 \\
    * - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 & \cos \theta_1 \sin \theta_2 & - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-212 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.point(0,1) = stheta1*stheta2;
    dcm.point(0,2) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.point(1,0) = stheta2*stheta3;
    dcm.point(1,1) = ctheta2;
    dcm.point(1,2) = -stheta2*ctheta3;
    dcm.point(2,0) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.point(2,1) = ctheta1*stheta2;
    dcm.point(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

    return dcm;
}


SO SO::from_eulerangles_body232(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \cos \theta_1 \sin \theta_2 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 \\
    * \sin \theta_2 \cos \theta_3 & \cos \theta_2 & \sin \theta_2 \sin \theta_3 \\
    * - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_2 & - \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-232 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.point(0,1) = -ctheta1*stheta2;
    dcm.point(0,2) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.point(1,0) = stheta2*ctheta3;
    dcm.point(1,1) = ctheta2;
    dcm.point(1,2) = stheta2*stheta3;
    dcm.point(2,0) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.point(2,1) = stheta1*stheta2;
    dcm.point(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

    return dcm;
}


SO SO::from_eulerangles_body313(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_2 \\
    * \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 - \cos \theta_1 \sin \theta_2 \\
    * \sin \theta_2 \sin \theta_3 & \sin \theta_2 \cos \theta_3 & \cos \theta_2
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-313 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.point(0,1) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.point(0,2) = stheta1*stheta2;
    dcm.point(1,0) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.point(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.point(1,2) = -ctheta1*stheta2;
    dcm.point(2,0) = stheta2*stheta3;
    dcm.point(2,1) = stheta2*ctheta3;
    dcm.point(2,2) = ctheta2;

    return dcm;
}


SO SO::from_eulerangles_body323(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 & \cos \theta_1 \sin \theta_2 \\
    * \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_1 \sin \theta_2 \\
    * - \sin \theta_2 \cos \theta_3 & \sin \theta_2 \sin \theta_3 & \cos \theta_2
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Body-323 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.point(0,1) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.point(0,2) = ctheta1*stheta2;
    dcm.point(1,0) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.point(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.point(1,2) = stheta1*stheta2;
    dcm.point(2,0) = -stheta2*ctheta3;
    dcm.point(2,1) = stheta2*stheta3;
    dcm.point(2,2) = ctheta2;

    return dcm;
}


SO SO::from_eulerangles_space123(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_2 \cos \theta_3 & - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 \\
    * \cos \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 \\
    * - \sin \theta_2 & \sin \theta_1 \cos \theta_2 & \cos \theta_1 \cos \theta_2
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-123 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta2*ctheta3;
    dcm.point(0,1) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(0,2) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.point(1,0) = ctheta2*stheta3;
    dcm.point(1,1) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
    dcm.point(1,2) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(2,0) = -stheta2;
    dcm.point(2,1) = stheta1*ctheta2;
    dcm.point(2,2) = ctheta1*ctheta2;

    return dcm;
}


SO SO::from_eulerangles_space231(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta 2 & - \sin \theta_2 & \sin \theta_1 \cos \theta_2 \\
    * \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_2 \cos \theta_3 & - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 \\
    * - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-231 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta2;
    dcm.point(0,1) = -stheta2;
    dcm.point(0,2) = stheta1*ctheta2;
    dcm.point(1,0) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.point(1,1) = ctheta2*ctheta3;
    dcm.point(1,2) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(2,0) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(2,1) = ctheta2*stheta3;
    dcm.point(2,2) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;

    return dcm;
}


SO SO::from_eulerangles_space312(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_3 + \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_2 \sin \theta_3 \\
    * \sin \theta_1 \cos \theta_2 & \cos \theta_1 \cos \theta_2 & - \sin \theta_2 \\
    * - \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_3 + \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_2 \cos \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-312 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
    dcm.point(0,1) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(0,2) = ctheta2*stheta3;
    dcm.point(1,0) = stheta1*ctheta2;
    dcm.point(1,1) = ctheta1*ctheta2;
    dcm.point(1,2) = -stheta2;
    dcm.point(2,0) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(2,1) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.point(2,2) = ctheta2*ctheta3;

    return dcm;
}


SO SO::from_eulerangles_space132(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 \\
    * \sin \theta_2 & \cos \theta_1 \cos \theta_2 & \sin \theta_1 \cos \theta_2 \\
    * - \cos \theta_2 \sin \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-132 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta2*ctheta3;
    dcm.point(0,1) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.point(0,2) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(1,0) = stheta2;
    dcm.point(1,1) = ctheta1*ctheta2;
    dcm.point(1,2) = -stheta1*ctheta2;
    dcm.point(2,0) = -ctheta2*stheta3;
    dcm.point(2,1) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(2,2) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;

    return dcm;
}


SO SO::from_eulerangles_space213(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \cos \theta_2 \sin \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 \\
    * \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_2 \cos \theta_3 & \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 \\
    * - \sin \theta_1 \cos \theta_2 & \sin \theta_2 & \cos \theta_1 \cos \theta_2
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-213 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
    dcm.point(0,1) = -ctheta2*stheta3;
    dcm.point(0,2) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(1,0) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(1,1) = ctheta2*ctheta3;
    dcm.point(1,2) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.point(2,0) = -stheta1*ctheta2;
    dcm.point(2,1) = stheta2;
    dcm.point(2,2) = ctheta1*ctheta2;

    return dcm;
}


SO SO::from_eulerangles_space321(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_2 & - \sin \theta_1 \cos \theta_2 & \sin \theta_2 \\
    * \sin \theta_1 \cos \theta_3 + \cos \theta_1 \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \sin \theta_2 \sin \theta_3 & - \cos \theta_2 \sin \theta_3 \\
    * \sin \theta_1 \sin \theta_3 - \cos \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \sin \theta_2 \cos \theta_3 & \cos \theta_2 \cos \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-321 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta2;
    dcm.point(0,1) = -stheta1*ctheta2;
    dcm.point(0,2) = stheta2;
    dcm.point(1,0) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.point(1,1) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
    dcm.point(1,2) = -ctheta2*stheta3;
    dcm.point(2,0) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.point(2,1) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.point(2,2) = ctheta2*ctheta3;

    return dcm;
}


SO SO::from_eulerangles_space121(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_2 & \sin \theta_1 \sin \theta_2 & \cos \theta_1 \sin \theta_2 \\
    * \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 \\
    * - \sin \theta_2 \cos \theta_3 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-121 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta2;
    dcm.point(0,1) = stheta1*stheta2;
    dcm.point(0,2) = ctheta1*stheta2;
    dcm.point(1,0) = stheta2*stheta3;
    dcm.point(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.point(1,2) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.point(2,0) = -stheta2*ctheta3;
    dcm.point(2,1) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.point(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

    return dcm;
}


SO SO::from_eulerangles_space131(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_2 & - \cos \theta_1 \sin \theta_2 & \sin \theta_1 \sin \theta_2 \\
    * \sin \theta_2 \cos \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 \\
    * \sin \theta_2 \sin \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-131 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta2;
    dcm.point(0,1) = -ctheta1*stheta2;
    dcm.point(0,2) = stheta1*stheta2;
    dcm.point(1,0) = stheta2*ctheta3;
    dcm.point(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.point(1,2) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.point(2,0) = stheta2*stheta3;
    dcm.point(2,1) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.point(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

    return dcm;
}


SO SO::from_eulerangles_space212(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_2 \sin \theta_3 & \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 \\
    * \sin \theta_1 \sin \theta_2 & \cos \theta_2 & - \cos \theta_1 \sin \theta_2 \\
    * - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 & \sin \theta_2 \cos \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-212 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.point(0,1) = stheta2*stheta3;
    dcm.point(0,2) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.point(1,0) = stheta1*stheta2;
    dcm.point(1,1) = ctheta2;
    dcm.point(1,2) = -ctheta1*stheta2;
    dcm.point(2,0) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.point(2,1) = stheta2*ctheta3;
    dcm.point(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

    return dcm;
}


SO SO::from_eulerangles_space232(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_2 \cos \theta_3 & \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 \\
    * \cos \theta_1 \sin \theta_2 & \cos \theta_2 & \sin \theta_1 \sin \theta_2 \\
    * - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-232 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.point(0,1) = -stheta2*ctheta3;
    dcm.point(0,2) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.point(1,0) = ctheta1*stheta2;
    dcm.point(1,1) = ctheta2;
    dcm.point(1,2) = stheta1*stheta2;
    dcm.point(2,0) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.point(2,1) = stheta2*stheta3;
    dcm.point(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

    return dcm;
}


SO SO::from_eulerangles_space313(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & - \sin \theta_1 \cos \theta_3 - \cos \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_2 \sin \theta_3 \\
    * \cos \theta_1 \sin \theta_3 + \sin \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \sin \theta_2 \cos \theta_3 \\
    * \sin \theta_1 \sin \theta_2 & \cos \theta_1 \sin \theta_2 & \cos \theta_2
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-313 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.point(0,1) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.point(0,2) = stheta2*stheta3;
    dcm.point(1,0) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.point(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.point(1,2) = -stheta2*ctheta3;
    dcm.point(2,0) = stheta1*stheta2;
    dcm.point(2,1) = ctheta1*stheta2;
    dcm.point(2,2) = ctheta2;

    return dcm;
}


SO SO::from_eulerangles_space323(const double theta1, const double theta2, const double theta3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * - \sin \theta_1 \sin \theta_3 + \cos \theta_1 \cos \theta_2 \cos \theta_3 & - \cos \theta_1 \sin \theta_3 - \sin \theta_1 \cos \theta_2 \cos \theta_3 & \sin \theta_2 \cos \theta_3 \\
    * \sin \theta_1 \cos \theta_3 + \cos \theta_1 \cos \theta_2 \sin \theta_3 & \cos \theta_1 \cos \theta_3 - \sin \theta_1 \cos \theta_2 \sin \theta_3 & \sin \theta_2 \sin \theta_3 \\
    * - \cos \theta_1 \sin \theta_2 & \sin \theta_1 \sin \theta_2 & \cos \theta_2
    * \end{bmatrix}\f}
    *
    * Transforms an Euler Angle Space-323 rotation sequence into a direction co-sine matrix.
    * @param[in] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[in] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[in] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    * @param[out] dcm The direction co-sine matrix, \f$dcm\f$.
    */

    SO dcm(3);

    const double stheta1 = std::sin(theta1);
    const double ctheta1 = std::cos(theta1);
    const double stheta2 = std::sin(theta2);
    const double ctheta2 = std::cos(theta2);
    const double stheta3 = std::sin(theta3);
    const double ctheta3 = std::cos(theta3);

    dcm.point(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.point(0,1) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.point(0,2) = stheta2*ctheta3;
    dcm.point(1,0) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.point(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.point(1,2) = stheta2*stheta3;
    dcm.point(2,0) = -ctheta1*stheta2;
    dcm.point(2,1) = stheta1*stheta2;
    dcm.point(2,2) = ctheta2;

    return dcm;
}

SO SO::from_quaternion(const double e0, const double e1, const double e2, const double e3)
{
    /*! \f{equation*}{ (\mathbb{R}^4) \rightarrow SO \f}
    *
    * Constructor instantiating a Quaternion as an \f$SO\f$ object.
    * 
    * Enables instatiation like:
    * 
    *     Lielab::domain::SO DCM0 = Lielab::domain::SO::from_quaternion(1.0, 0.0, 0.0, 0.0);
    *
    *
    * \f{equation*}{dcm = \begin{bmatrix}
    * e_0^2 + e_1^2 - e_2^2 - e_3^2 & 2(e1 e2 - e3 e0) & 2(e1 e3 + e2 e0) \\
    * 2(e1 e2 + e3 e0) & e_0^2 - e_1^2 + e_2^2 - e_3^2 & 2(e2 e3 - e1 e0) \\
    * 2(e1 e3 - e2 e0) & 2(e2 e3 + e1 e0) & e_0^2 - e_1^2 - e_2^2 + e_3^2
    * \end{bmatrix}\f}
    * 
    * @param[out] dcm An SO object representing the input quaternion.
    *
    * TODO: Citation needed
    *
    */


    SO dcm(3);

    dcm.point(0,0) = e0*e0 + e1*e1 - e2*e2 - e3*e3;
    dcm.point(0,1) = 2.0*e1*e2 - 2.0*e3*e0;
    dcm.point(0,2) = 2.0*e1*e3 + 2.0*e2*e0;
    dcm.point(1,0) = 2.0*e1*e2 + 2.0*e3*e0;
    dcm.point(1,1) = e0*e0 - e1*e1 + e2*e2 - e3*e3;
    dcm.point(1,2) = 2.0*e2*e3 - 2.0*e1*e0;
    dcm.point(2,0) = 2.0*e1*e3 - 2.0*e2*e0;
    dcm.point(2,1) = 2.0*e2*e3 + 2.0*e1*e0;
    dcm.point(2,2) = e0*e0 - e1*e1 - e2*e2 + e3*e3;

    return dcm;
}


SO SO::from_rodriguesvector(const double g1, const double g2, const double g3)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathbb{R}, \mathbb{R}) \rightarrow SO(3) \f}
    *
    * \f{equation*}{dcm = \frac{1}{1 + g_1^2 + g_2^2 + g_3^2}\begin{bmatrix}
    * 1 + g_1^2 - g_2^2 - g_3^2 & 2(g_1 g_2 + g_3) & 2(g_1 g_3 - g_2) \\
    * 2(g_1 g_2 - g_3) & 1 - g_1^2 + g_2^2 - g_3^2 & 2(g_2 g_3 + g_1) \\
    * 2(g_1 g_3 + g_2) & 2(g_2 g_3 - g_1) & 1 - g_1^2 - g_2^2 + g_3^2
    * \end{bmatrix}\f}
    * 
    * Transforms a Rodrigues (Gibbs) vector into a direction co-sine matrix.
    * @param[in] g1 The 1st element in Rodrigues representation.
    * @param[in] g2 The 2nd element in Rodrigues representation.
    * @param[in] g3 The 3rd element in Rodrigues representation.
    * @param[out] dcm The direction co-sine matrix.
    * 
    * TODO: Citation needed
    * TODO: Check conventions. This might be transposed from what we use elsewhere.
    * 
    */

    SO dcm(3);

    const double g12 = std::pow(g1, 2.0);
    const double g22 = std::pow(g2, 2.0);
    const double g32 = std::pow(g3, 2.0);
    const double mul = 1.0/(1.0 + g12 + g22 + g32);

    dcm.point(0,0) = mul*(1.0 + g12 - g22 - g32);
    dcm.point(0,1) = mul*(2.0*(g1*g2 + g3));
    dcm.point(0,2) = mul*(2.0*(g1*g3 - g2));
    dcm.point(1,0) = mul*(2.0*(g1*g2 - g3));
    dcm.point(1,1) = mul*(1.0 - g12 + g22 - g32);
    dcm.point(1,2) = mul*(2.0*(g2*g3 + g1));
    dcm.point(2,0) = mul*(2.0*(g1*g3 + g2));
    dcm.point(2,1) = mul*(2.0*(g2*g3 - g1));
    dcm.point(2,2) = mul*(1.0 - g12 - g22 + g32);

    return dcm;
}

SO SO::from_SU2(const SU& other)
{
    /*!
        * Transforms an SU(2) object to a Direction Cosine object.
        *
        * @param[in] q A quaternion as an SU object of shape 2.
        * @param[out] dcm A direction cosine as an SO object of shape 3.
        */

    if (other.get_shape() != 2)
    {
        throw std::domain_error("SO::from_SU2: Expected input shape 2. Got " + std::to_string(other.get_shape()) + ".");
    }

    const auto [q0, q1, q2, q3] = other.to_quaternion();

    SO dcm(3);
    dcm.point(0,0) = 1.0 - 2.0*q2*q2 - 2.0*q3*q3;
    dcm.point(0,1) = 2.0*q1*q2 - 2.0*q3*q0;
    dcm.point(0,2) = 2.0*q1*q3 + 2.0*q2*q0;
    dcm.point(1,0) = 2.0*q1*q2 + 2.0*q3*q0;
    dcm.point(1,1) = 1.0 - 2.0*q1*q1 - 2.0*q3*q3;
    dcm.point(1,2) = 2.0*q2*q3 - 2.0*q1*q0;
    dcm.point(2,0) = 2.0*q1*q3 - 2.0*q2*q0;
    dcm.point(2,1) = 2.0*q2*q3 + 2.0*q1*q0;
    dcm.point(2,2) = 1.0 - 2.0*q1*q1 - 2.0*q2*q2;

    return dcm;
}

std::string SO::to_string() const
{
    return "SO(" + std::to_string(this->get_shape()) + ")";
}

int SO::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    return this->get_shape() * (this->get_shape() - 1) / 2;
}

int SO::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
    *
    * Gets the size of the data representation.
    */

    if (this->get_shape() <= 1) return 0;

    return static_cast<int>(std::pow(this->get_shape(), 2));
}

bool SO::is_abelian() const
{
    return false;
}

int SO::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return static_cast<int>(this->point.rows());
}

SO::point_t SO::get_point() const
{
    return this->point;
}

Eigen::VectorXd SO::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */

    if (this->get_shape() <= 1) return Eigen::VectorXd::Zero(0);

    return this->point.reshaped<Eigen::RowMajor>();
}

void SO::unserialize(const Eigen::VectorXd& serialized)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the GL object from a serialized vector.
    */

    const int vdim = static_cast<int>(serialized.size());
    const int max_ind = std::min(this->get_size(), vdim);

    for (int vind = 0; vind < max_ind; vind++)
    {
        const int row = vind / this->get_shape();
        const int col = vind % this->get_shape();
        this->point(row, col) = serialized(vind);
    }
}

void SO::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

SO::matrix_t SO::get_matrix() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times n} \f}
    * 
    * Returns a matrix representation.
    * 
    * Formerly called "get_ados_representation()".
    * 
    * Ado, Igor D. "Note on the representation of finite continuous groups by
    *               means of linear substitutions, Izv. Fiz." Mat. Obsch.(Kazan)
    *               7.1 (1935): 935.
    * 
    * Ado, Igor D. "The representation of Lie algebras by matrices." Uspekhi
    *               Matematicheskikh Nauk 2.6 (1947): 159-173.
    */

    return this->point;
}

SO::field_t SO::operator()(const int index1, const int index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int shape = this->get_shape();
    if (shape == 0) return nan;

    // If input index is negative, index from the back of the array
    const int _index1 = (index1 < 0) ? shape + index1 : index1;
    const int _index2 = (index2 < 0) ? shape + index2 : index2;

    // Error check for out of bounds
    if (_index1 < 0) return nan;
    if (_index1 >= shape) return nan;
    if (_index2 < 0) return nan;
    if (_index2 >= shape) return nan;

    return this->point(_index1, _index2);
}

std::array<double, 3> SO::to_eulerangles_body123() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(this->point(0,2));

    if (std::abs(this->point(0,2) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(1,0);
        const double ctheta1 = -this->point(2,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(0,2) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(1,0);
        const double ctheta1 = this->point(2,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = -this->point(1,2) / ctheta2;
    const double ctheta1 = this->point(2,2) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = -this->point(0,1) / ctheta2;
    const double ctheta3 = this->point(0,0) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body231() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(this->point(1,0));

    if (std::abs(this->point(1,0) - 1) < 1e-14)
    {
        const double stheta1 = this->point(2,1);
        const double ctheta1 = -this->point(0,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(1,0) + 1) < 1e-14)
    {
        const double stheta1 = -this->point(2,1);
        const double ctheta1 = this->point(0,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = -this->point(2,0) / ctheta2;
    const double ctheta1 = this->point(0,0) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = -this->point(1,2) / ctheta2;
    const double ctheta3 = this->point(1,1) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body312() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(this->point(2,1));

    if (std::abs(this->point(2,1) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(0,2);
        const double ctheta1 = -this->point(1,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(2,1) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(0,2);
        const double ctheta1 = this->point(1,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = -this->point(0,1) / ctheta2;
    const double ctheta1 = this->point(1,1) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = -this->point(2,0) / ctheta2;
    const double ctheta3 = this->point(2,2) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body132() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(-this->point(0,1));

    if (std::abs(this->point(0,1) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(2,0);
        const double ctheta1 = -this->point(1,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(0,1) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(2,0);
        const double ctheta1 = this->point(1,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = this->point(2,1) / ctheta2;
    const double ctheta1 = this->point(1,1) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(0,2) / ctheta2;
    const double ctheta3 = this->point(0,0) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body213() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(-this->point(1,2));

    if (std::abs(this->point(1,2) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(0,1);
        const double ctheta1 = -this->point(2,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(1,2) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(0,1);
        const double ctheta1 = this->point(2,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = this->point(0,2) / ctheta2;
    const double ctheta1 = this->point(2,2) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(1,0) / ctheta2;
    const double ctheta3 = this->point(1,1) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body321() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(-this->point(2,0));

    if (std::abs(this->point(2,0) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(1,2);
        const double ctheta1 = -this->point(0,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(2,0) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(1,2);
        const double ctheta1 = this->point(0,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = this->point(1,0) / ctheta2;
    const double ctheta1 = this->point(0,0) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(2,1) / ctheta2;
    const double ctheta3 = this->point(2,2) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body121() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(0,0));

    if (std::abs(this->point(0,0) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(1,2);
        const double ctheta1 = this->point(2,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(0,0) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(1,2);
        const double ctheta1 = -this->point(2,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(1,0) / stheta2;
    const double ctheta1 = -this->point(2,0) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(0,1) / stheta2;
    const double ctheta3 = this->point(0,2) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body131() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(0,0));

    if (std::abs(this->point(0,0) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(2,1);
        const double ctheta1 = this->point(1,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(0,0) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(2,1);
        const double ctheta1 = -this->point(1,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(2,0) / stheta2;
    const double ctheta1 = this->point(1,0) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(0,2) / stheta2;
    const double ctheta3 = -this->point(0,1) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body212() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(1,1));

    if (std::abs(this->point(1,1) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(0,2);
        const double ctheta1 = this->point(2,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(1,1) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(0,2);
        const double ctheta1 = -this->point(2,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(0,1) / stheta2;
    const double ctheta1 = this->point(2,1) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(1,0) / stheta2;
    const double ctheta3 = -this->point(1,2) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body232() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(1,1));

    if (std::abs(this->point(1,1) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(2,0);
        const double ctheta1 = this->point(0,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(1,1) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(2,0);
        const double ctheta1 = -this->point(0,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(2,1) / stheta2;
    const double ctheta1 = -this->point(0,1) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(1,2) / stheta2;
    const double ctheta3 = this->point(1,0) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body313() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(2,2));

    if (std::abs(this->point(2,2) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(0,1);
        const double ctheta1 = this->point(1,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(2,2) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(0,1);
        const double ctheta1 = -this->point(1,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(0,2) / stheta2;
    const double ctheta1 = -this->point(1,2) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(2,0) / stheta2;
    const double ctheta3 = this->point(2,1) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_body323() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(2,2));

    if (std::abs(this->point(2,2) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(1,0);
        const double ctheta1 = this->point(0,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(2,2) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(1,0);
        const double ctheta1 = -this->point(0,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(1,2) / stheta2;
    const double ctheta1 = this->point(0,2) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(2,1) / stheta2;
    const double ctheta3 = -this->point(2,0) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space123() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(-this->point(2,0));

    if (std::abs(this->point(2,0) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(0,1);
        const double ctheta1 = -this->point(0,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(2,0) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(0,1);
        const double ctheta1 = this->point(0,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = this->point(2,1) / ctheta2;
    const double ctheta1 = this->point(2,2) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(1,0) / ctheta2;
    const double ctheta3 = this->point(0,0) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space231() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(-this->point(0,1));

    if (std::abs(this->point(0,1) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(1,2);
        const double ctheta1 = -this->point(1,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(0,1) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(1,2);
        const double ctheta1 = this->point(1,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = this->point(0,2) / ctheta2;
    const double ctheta1 = this->point(0,0) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(2,1) / ctheta2;
    const double ctheta3 = this->point(1,1) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space312() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(-this->point(1,2));

    if (std::abs(this->point(1,2) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(2,0);
        const double ctheta1 = -this->point(2,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(1,2) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(2,0);
        const double ctheta1 = this->point(2,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = this->point(1,0) / ctheta2;
    const double ctheta1 = this->point(1,1) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(0,2) / ctheta2;
    const double ctheta3 = this->point(2,2) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space132() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(this->point(1,0));

    if (std::abs(this->point(1,0) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(0,2);
        const double ctheta1 = -this->point(0,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(1,0) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(0,2);
        const double ctheta1 = this->point(0,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = -this->point(1,2) / ctheta2;
    const double ctheta1 = this->point(1,1) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = -this->point(2,0) / ctheta2;
    const double ctheta3 = this->point(0,0) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}



std::array<double, 3> SO::to_eulerangles_space213() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(this->point(2,1));

    if (std::abs(this->point(2,1) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(1,0);
        const double ctheta1 = -this->point(1,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(2,1) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(1,0);
        const double ctheta1 = this->point(1,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = -this->point(2,0) / ctheta2;
    const double ctheta1 = this->point(2,2) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = -this->point(0,1) / ctheta2;
    const double ctheta3 = this->point(1,1) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space321() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::asin(this->point(0,2));

    if (std::abs(this->point(0,2) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(2,1);
        const double ctheta1 = -this->point(2,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(0,2) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(2,1);
        const double ctheta1 = this->point(2,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double ctheta2 = std::cos(theta2);

    const double stheta1 = -this->point(0,1) / ctheta2;
    const double ctheta1 = this->point(0,0) / ctheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = -this->point(1,2) / ctheta2;
    const double ctheta3 = this->point(2,2) / ctheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space121() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(0,0));

    if (std::abs(this->point(0,0) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(2,1);
        const double ctheta1 = this->point(2,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(0,0) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(2,1);
        const double ctheta1 = -this->point(2,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(0,1) / stheta2;
    const double ctheta1 = this->point(0,2) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(1,0) / stheta2;
    const double ctheta3 = -this->point(2,0) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space131() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(0,0));

    if (std::abs(this->point(0,0) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(1,2);
        const double ctheta1 = this->point(1,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(0,0) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(1,2);
        const double ctheta1 = -this->point(1,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(0,2) / stheta2;
    const double ctheta1 = -this->point(0,1) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(2,0) / stheta2;
    const double ctheta3 = this->point(1,0) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space212() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(1,1));

    if (std::abs(this->point(1,1) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(2,0);
        const double ctheta1 = this->point(2,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(1,1) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(2,0);
        const double ctheta1 = -this->point(2,2);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(1,0) / stheta2;
    const double ctheta1 = -this->point(1,2) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(0,1) / stheta2;
    const double ctheta3 = this->point(2,1) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space232() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(1,1));

    if (std::abs(this->point(1,1) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(0,2);
        const double ctheta1 = this->point(0,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(1,1) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(0,2);
        const double ctheta1 = -this->point(0,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(1,2) / stheta2;
    const double ctheta1 = this->point(1,0) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(2,1) / stheta2;
    const double ctheta3 = -this->point(0,1) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space313() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(2,2));

    if (std::abs(this->point(2,2) - 1.0) < 1e-14)
    {
        const double stheta1 = this->point(1,0);
        const double ctheta1 = this->point(1,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(2,2) + 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(1,0);
        const double ctheta1 = -this->point(1,1);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(2,0) / stheta2;
    const double ctheta1 = this->point(2,1) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(0,2) / stheta2;
    const double ctheta3 = -this->point(1,2) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


std::array<double, 3> SO::to_eulerangles_space323() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->point The direction co-sine matrix, \f$this->point\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double theta2 = std::acos(this->point(2,2));

    if (std::abs(this->point(2,2) - 1.0) < 1e-14)
    {
        const double stheta1 = -this->point(0,1);
        const double ctheta1 = this->point(0,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->point(2,2) + 1.0) < 1e-14)
    {
        const double stheta1 = this->point(0,1);
        const double ctheta1 = -this->point(0,0);
        const double theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const double stheta2 = std::sin(theta2);

    const double stheta1 = this->point(2,1) / stheta2;
    const double ctheta1 = -this->point(2,0) / stheta2;
    const double theta1 = std::atan2(stheta1, ctheta1);

    const double stheta3 = this->point(1,2) / stheta2;
    const double ctheta3 = this->point(0,2) / stheta2;
    const double theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

std::array<double, 4> SO::to_quaternion() const
{
    /*! \f{equation*}{ (SO) \rightarrow \mathbb{R}^4 \f}
    *
    * tbd...
    *
    * TODO: Citation needed
    *
    */

    lielab_assert(this->get_shape() == 3, "Expected shape 3. Got " + std::to_string(this->get_shape()) + ".");

    const double e[] = {0.5*std::sqrt(1.0 + this->point(0,0) + this->point(1,1) + this->point(2,2)),
                        0.5*std::sqrt(1.0 + this->point(0,0) - this->point(1,1) - this->point(2,2)),
                        0.5*std::sqrt(1.0 - this->point(0,0) + this->point(1,1) - this->point(2,2)),
                        0.5*std::sqrt(1.0 - this->point(0,0) - this->point(1,1) + this->point(2,2))};
    
    int ind = 0;
    double max = 0;

    for (int ii = 0; ii < 4; ii++)
    {
        if (max < e[ii])
        {
            ind = ii;
            max = e[ii];
        }
    }

    if (ind == 0)
    {
        return {e[0], (this->point(2,1) - this->point(1,2))/(4.0*e[0]), (this->point(0,2) - this->point(2,0))/(4.0*e[0]), (this->point(1,0) - this->point(0,1))/(4.0*e[0])};
    }
    else if (ind == 1)
    {
        return {(this->point(2,1) - this->point(1,2))/(4.0*e[1]), e[1], (this->point(1,0) + this->point(0,1))/(4.0*e[1]), (this->point(0,2) + this->point(2,0))/(4.0*e[1])};
    }
    else if (ind == 2)
    {
        return {(this->point(0,2) - this->point(2,0))/(4.0*e[2]), (this->point(1,0) + this->point(0,1))/(4.0*e[2]), e[2], (this->point(2,1) + this->point(1,2))/(4.0*e[2])};
    }

    // ind == 3 should be the only other option by this point.
    assert(ind == 3);

    return {(this->point(1,0) - this->point(0,1))/(4.0*e[3]), (this->point(0,2) + this->point(2,0))/(4.0*e[3]), (this->point(2,1) + this->point(1,2))/(4.0*e[3]), e[3]};
}

std::array<double, 3> SO::to_gibbs() const
{
    /*!
    * Transforms a direction co-sine matrics into a Gibbs (Rodrigues) vector.
    * @param[in] dcm The direction co-sine matrix.
    * @param[out] gibbs The rotation in Gibbs representation.
    */

    const double mul = 1.0/(1.0 + this->point(0,0) + this->point(1,1) + this->point(2,2));
    
    const double g0 = mul*(this->point(1,2) - this->point(2,1));
    const double g1 = mul*(this->point(2,0) - this->point(0,2));
    const double g2 = mul*(this->point(0,1) - this->point(1,0));

    return {g0, g1, g2};
}

SO SO::operator*(const SO& other) const
{
    /*! \f{equation*}{ (SO, SO) \rightarrow SO \f}
    *
    * Group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return SO(this->point * other.point);
}

SO& SO::operator*=(const SO& other)
{
    /*! \f{equation*}{ (SO, SO) \rightarrow SO \f}
    *
    * In place group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point *= other.point;
    return *this;
}

SO SO::inverse() const
{
    /*! \f{equation*}{ (SO) \rightarrow SO \f}
    * 
    * Returns the inverse.
    */

    return SO(this->point.transpose());
}

so SO::project_onto_tangent_space(const Eigen::MatrixXd& vector)
{
    // TODO: Error check here
    const Eigen::MatrixXd o1 = this->inverse().point*vector;
    const Eigen::MatrixXd m = 0.5*(o1 - o1.transpose());
    return so(m);
}

}
