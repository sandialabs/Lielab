#ifndef LIELAB_DOMAIN_SO_TPP
#define LIELAB_DOMAIN_SO_TPP

#include "SO.hpp"

#include "SU.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <cmath>

namespace Lielab::domain
{

template<typename OtherDerived>
SO::SO(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{equation}{ (\mathbb{R}^{n \times n}) \rightarrow SO \f}
    *
    * Constructor instantiating an \f$SO\f$ object from an
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
SO & SO::operator=(const Eigen::MatrixBase<OtherDerived> & other)
{
    /*! \f{equation}{ SO := \mathbb{R}^{n \times n} \f}
    * 
    * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `SO`.
    */

    if (other.rows() != other.cols())
    {
        throw Lielab::utils::Error("Size of the matrix must be square.");
    }
    
    this->data = Eigen::MatrixXd(other);
    this->_shape = this->data.rows();
    return *this;
}

/*
* Additional static initializers. Not a part of the core Lie group, but are convenient.
*/
template <typename T>
SO SO::from_eulerangles_body123(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta2*ctheta3;
    dcm.data(0,1) = -ctheta2*stheta3;
    dcm.data(0,2) = stheta2;
    dcm.data(1,0) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(1,1) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
    dcm.data(1,2) = -stheta1*ctheta2;
    dcm.data(2,0) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.data(2,1) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(2,2) = ctheta1*ctheta2;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body231(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta2;
    dcm.data(0,1) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.data(0,2) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(1,0) = stheta2;
    dcm.data(1,1) = ctheta2*ctheta3;
    dcm.data(1,2) = -ctheta2*stheta3;
    dcm.data(2,0) = -stheta1*ctheta2;
    dcm.data(2,1) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(2,2) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body312(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
    dcm.data(0,1) = -stheta1*ctheta2;
    dcm.data(0,2) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(1,0) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(1,1) = ctheta1*ctheta2;
    dcm.data(1,2) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.data(2,0) = -ctheta2*stheta3;
    dcm.data(2,1) = stheta2;
    dcm.data(2,2) = ctheta2*ctheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body132(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta2*ctheta3;
    dcm.data(0,1) = -stheta2;
    dcm.data(0,2) = ctheta2*stheta3;
    dcm.data(1,0) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.data(1,1) = ctheta1*ctheta2;
    dcm.data(1,2) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(2,0) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(2,1) = stheta1*ctheta2;
    dcm.data(2,2) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body213(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
    dcm.data(0,1) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(0,2) = stheta1*ctheta2;
    dcm.data(1,0) = ctheta2*stheta3;
    dcm.data(1,1) = ctheta2*ctheta3;
    dcm.data(1,2) = -stheta2;
    dcm.data(2,0) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(2,1) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.data(2,2) = ctheta1*ctheta2;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body321(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta2;
    dcm.data(0,1) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(0,2) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.data(1,0) = stheta1*ctheta2;
    dcm.data(1,1) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
    dcm.data(1,2) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(2,0) = -stheta2;
    dcm.data(2,1) = ctheta2*stheta3;
    dcm.data(2,2) = ctheta2*ctheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body121(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta2;
    dcm.data(0,1) = stheta2*stheta3;
    dcm.data(0,2) = stheta2*ctheta3;
    dcm.data(1,0) = stheta1*stheta2;
    dcm.data(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.data(1,2) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.data(2,0) = -ctheta1*stheta2;
    dcm.data(2,1) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.data(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body131(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta2;
    dcm.data(0,1) = -stheta2*ctheta3;
    dcm.data(0,2) = stheta2*stheta3;
    dcm.data(1,0) = ctheta1*stheta2;
    dcm.data(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.data(1,2) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.data(2,0) = stheta1*stheta2;
    dcm.data(2,1) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.data(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body212(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.data(0,1) = stheta1*stheta2;
    dcm.data(0,2) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.data(1,0) = stheta2*stheta3;
    dcm.data(1,1) = ctheta2;
    dcm.data(1,2) = -stheta2*ctheta3;
    dcm.data(2,0) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.data(2,1) = ctheta1*stheta2;
    dcm.data(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body232(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.data(0,1) = -ctheta1*stheta2;
    dcm.data(0,2) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.data(1,0) = stheta2*ctheta3;
    dcm.data(1,1) = ctheta2;
    dcm.data(1,2) = stheta2*stheta3;
    dcm.data(2,0) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.data(2,1) = stheta1*stheta2;
    dcm.data(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body313(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.data(0,1) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.data(0,2) = stheta1*stheta2;
    dcm.data(1,0) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.data(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.data(1,2) = -ctheta1*stheta2;
    dcm.data(2,0) = stheta2*stheta3;
    dcm.data(2,1) = stheta2*ctheta3;
    dcm.data(2,2) = ctheta2;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_body323(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.data(0,1) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.data(0,2) = ctheta1*stheta2;
    dcm.data(1,0) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.data(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.data(1,2) = stheta1*stheta2;
    dcm.data(2,0) = -stheta2*ctheta3;
    dcm.data(2,1) = stheta2*stheta3;
    dcm.data(2,2) = ctheta2;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space123(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta2*ctheta3;
    dcm.data(0,1) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(0,2) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.data(1,0) = ctheta2*stheta3;
    dcm.data(1,1) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
    dcm.data(1,2) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(2,0) = -stheta2;
    dcm.data(2,1) = stheta1*ctheta2;
    dcm.data(2,2) = ctheta1*ctheta2;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space231(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta2;
    dcm.data(0,1) = -stheta2;
    dcm.data(0,2) = stheta1*ctheta2;
    dcm.data(1,0) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.data(1,1) = ctheta2*ctheta3;
    dcm.data(1,2) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(2,0) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(2,1) = ctheta2*stheta3;
    dcm.data(2,2) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space312(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta3 + stheta1*stheta2*stheta3;
    dcm.data(0,1) = -stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(0,2) = ctheta2*stheta3;
    dcm.data(1,0) = stheta1*ctheta2;
    dcm.data(1,1) = ctheta1*ctheta2;
    dcm.data(1,2) = -stheta2;
    dcm.data(2,0) = -ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(2,1) = stheta1*stheta3 + ctheta1*stheta2*ctheta3;
    dcm.data(2,2) = ctheta2*ctheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space132(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta2*ctheta3;
    dcm.data(0,1) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.data(0,2) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(1,0) = stheta2;
    dcm.data(1,1) = ctheta1*ctheta2;
    dcm.data(1,2) = -stheta1*ctheta2;
    dcm.data(2,0) = -ctheta2*stheta3;
    dcm.data(2,1) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(2,2) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space213(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
    dcm.data(0,1) = -ctheta2*stheta3;
    dcm.data(0,2) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(1,0) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(1,1) = ctheta2*ctheta3;
    dcm.data(1,2) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.data(2,0) = -stheta1*ctheta2;
    dcm.data(2,1) = stheta2;
    dcm.data(2,2) = ctheta1*ctheta2;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space321(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta2;
    dcm.data(0,1) = -stheta1*ctheta2;
    dcm.data(0,2) = stheta2;
    dcm.data(1,0) = stheta1*ctheta3 + ctheta1*stheta2*stheta3;
    dcm.data(1,1) = ctheta1*ctheta3 - stheta1*stheta2*stheta3;
    dcm.data(1,2) = -ctheta2*stheta3;
    dcm.data(2,0) = stheta1*stheta3 - ctheta1*stheta2*ctheta3;
    dcm.data(2,1) = ctheta1*stheta3 + stheta1*stheta2*ctheta3;
    dcm.data(2,2) = ctheta2*ctheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space121(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta2;
    dcm.data(0,1) = stheta1*stheta2;
    dcm.data(0,2) = ctheta1*stheta2;
    dcm.data(1,0) = stheta2*stheta3;
    dcm.data(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.data(1,2) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.data(2,0) = -stheta2*ctheta3;
    dcm.data(2,1) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.data(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space131(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta2;
    dcm.data(0,1) = -ctheta1*stheta2;
    dcm.data(0,2) = stheta1*stheta2;
    dcm.data(1,0) = stheta2*ctheta3;
    dcm.data(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.data(1,2) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.data(2,0) = stheta2*stheta3;
    dcm.data(2,1) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.data(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space212(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.data(0,1) = stheta2*stheta3;
    dcm.data(0,2) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.data(1,0) = stheta1*stheta2;
    dcm.data(1,1) = ctheta2;
    dcm.data(1,2) = -ctheta1*stheta2;
    dcm.data(2,0) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.data(2,1) = stheta2*ctheta3;
    dcm.data(2,2) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space232(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.data(0,1) = -stheta2*ctheta3;
    dcm.data(0,2) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.data(1,0) = ctheta1*stheta2;
    dcm.data(1,1) = ctheta2;
    dcm.data(1,2) = stheta1*stheta2;
    dcm.data(2,0) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.data(2,1) = stheta2*stheta3;
    dcm.data(2,2) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space313(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.data(0,1) = -stheta1*ctheta3 - ctheta1*ctheta2*stheta3;
    dcm.data(0,2) = stheta2*stheta3;
    dcm.data(1,0) = ctheta1*stheta3 + stheta1*ctheta2*ctheta3;
    dcm.data(1,1) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.data(1,2) = -stheta2*ctheta3;
    dcm.data(2,0) = stheta1*stheta2;
    dcm.data(2,1) = ctheta1*stheta2;
    dcm.data(2,2) = ctheta2;

    return dcm;
}

template <typename T>
SO SO::from_eulerangles_space323(const T theta1, const T theta2, const T theta3)
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

    const T stheta1 = std::sin(theta1);
    const T ctheta1 = std::cos(theta1);
    const T stheta2 = std::sin(theta2);
    const T ctheta2 = std::cos(theta2);
    const T stheta3 = std::sin(theta3);
    const T ctheta3 = std::cos(theta3);

    dcm.data(0,0) = -stheta1*stheta3 + ctheta1*ctheta2*ctheta3;
    dcm.data(0,1) = -ctheta1*stheta3 - stheta1*ctheta2*ctheta3;
    dcm.data(0,2) = stheta2*ctheta3;
    dcm.data(1,0) = stheta1*ctheta3 + ctheta1*ctheta2*stheta3;
    dcm.data(1,1) = ctheta1*ctheta3 - stheta1*ctheta2*stheta3;
    dcm.data(1,2) = stheta2*stheta3;
    dcm.data(2,0) = -ctheta1*stheta2;
    dcm.data(2,1) = stheta1*stheta2;
    dcm.data(2,2) = ctheta2;

    return dcm;
}

template <typename T>
SO SO::from_quaternion(const T e0, const T e1, const T e2, const T e3)
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

    dcm.data(0,0) = e0*e0 + e1*e1 - e2*e2 - e3*e3;
    dcm.data(0,1) = 2.0*e1*e2 - 2.0*e3*e0;
    dcm.data(0,2) = 2.0*e1*e3 + 2.0*e2*e0;
    dcm.data(1,0) = 2.0*e1*e2 + 2.0*e3*e0;
    dcm.data(1,1) = e0*e0 - e1*e1 + e2*e2 - e3*e3;
    dcm.data(1,2) = 2.0*e2*e3 - 2.0*e1*e0;
    dcm.data(2,0) = 2.0*e1*e3 - 2.0*e2*e0;
    dcm.data(2,1) = 2.0*e2*e3 + 2.0*e1*e0;
    dcm.data(2,2) = e0*e0 - e1*e1 - e2*e2 + e3*e3;

    return dcm;
}

template <typename T>
SO SO::from_rodriguesvector(const T g1, const T g2, const T g3)
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

    const T g12 = std::pow(g1, T(2.0));
    const T g22 = std::pow(g2, T(2.0));
    const T g32 = std::pow(g3, T(2.0));
    const T mul = 1.0/(1.0 + g12 + g22 + g32);

    dcm.data(0,0) = mul*(1.0 + g12 - g22 - g32);
    dcm.data(0,1) = mul*(2.0*(g1*g2 + g3));
    dcm.data(0,2) = mul*(2.0*(g1*g3 - g2));
    dcm.data(1,0) = mul*(2.0*(g1*g2 - g3));
    dcm.data(1,1) = mul*(1.0 - g12 + g22 - g32);
    dcm.data(1,2) = mul*(2.0*(g2*g3 + g1));
    dcm.data(2,0) = mul*(2.0*(g1*g3 + g2));
    dcm.data(2,1) = mul*(2.0*(g2*g3 - g1));
    dcm.data(2,2) = mul*(1.0 - g12 - g22 + g32);

    return dcm;
}

/*
* Additional outputs. Not a part of the core Lie group, but are convenient.
*/

template <typename T>
std::array<T, 4> SO::to_quaternion() const
{
    /*! \f{equation*}{ (SO) \rightarrow \mathbb{R}^4 \f}
    *
    * tbd...
    *
    * TODO: Citation needed
    *
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_quaternion: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const double e[] = {0.5*std::sqrt(1.0 + this->data(0,0) + this->data(1,1) + this->data(2,2)),
                        0.5*std::sqrt(1.0 + this->data(0,0) - this->data(1,1) - this->data(2,2)),
                        0.5*std::sqrt(1.0 - this->data(0,0) + this->data(1,1) - this->data(2,2)),
                        0.5*std::sqrt(1.0 - this->data(0,0) - this->data(1,1) + this->data(2,2))};
    
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
        return {e[0], (this->data(2,1) - this->data(1,2))/(4.0*e[0]), (this->data(0,2) - this->data(2,0))/(4.0*e[0]), (this->data(1,0) - this->data(0,1))/(4.0*e[0])};
    }
    else if (ind == 1)
    {
        return {(this->data(2,1) - this->data(1,2))/(4.0*e[1]), e[1], (this->data(1,0) + this->data(0,1))/(4.0*e[1]), (this->data(0,2) + this->data(2,0))/(4.0*e[1])};
    }
    else if (ind == 2)
    {
        return {(this->data(0,2) - this->data(2,0))/(4.0*e[2]), (this->data(1,0) + this->data(0,1))/(4.0*e[2]), e[2], (this->data(2,1) + this->data(1,2))/(4.0*e[2])};
    }

    // ind == 3 should be the only other option by this point.
    assert(ind == 3);

    return {(this->data(1,0) - this->data(0,1))/(4.0*e[3]), (this->data(0,2) + this->data(2,0))/(4.0*e[3]), (this->data(2,1) + this->data(1,2))/(4.0*e[3]), e[3]};
}

template <typename T>
std::array<T, 3> SO::to_gibbs() const
{
    /*!
    * Transforms a direction co-sine matrics into a Gibbs (Rodrigues) vector.
    * @param[in] dcm The direction co-sine matrix.
    * @param[out] gibbs The rotation in Gibbs representation.
    */

    const T mul = 1.0/(1.0 + this->data(0,0) + this->data(1,1) + this->data(2,2));
    
    const T g0 = mul*(this->data(1,2) - this->data(2,1));
    const T g1 = mul*(this->data(2,0) - this->data(0,2));
    const T g2 = mul*(this->data(0,1) - this->data(1,0));

    return {g0, g1, g2};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body123() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(this->data(0,2));

    if (std::abs(this->data(0,2) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(1,0);
        const T ctheta1 = -this->data(2,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(0,2) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(1,0);
        const T ctheta1 = this->data(2,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = -this->data(1,2) / ctheta2;
    const T ctheta1 = this->data(2,2) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = -this->data(0,1) / ctheta2;
    const T ctheta3 = this->data(0,0) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body231() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(this->data(1,0));

    if (std::abs(this->data(1,0) - 1) < 1e-14)
    {
        const T stheta1 = this->data(2,1);
        const T ctheta1 = -this->data(0,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(1,0) + 1) < 1e-14)
    {
        const T stheta1 = -this->data(2,1);
        const T ctheta1 = this->data(0,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = -this->data(2,0) / ctheta2;
    const T ctheta1 = this->data(0,0) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = -this->data(1,2) / ctheta2;
    const T ctheta3 = this->data(1,1) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body312() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(this->data(2,1));

    if (std::abs(this->data(2,1) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(0,2);
        const T ctheta1 = -this->data(1,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(2,1) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(0,2);
        const T ctheta1 = this->data(1,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = -this->data(0,1) / ctheta2;
    const T ctheta1 = this->data(1,1) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = -this->data(2,0) / ctheta2;
    const T ctheta3 = this->data(2,2) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body132() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(-this->data(0,1));

    if (std::abs(this->data(0,1) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(2,0);
        const T ctheta1 = -this->data(1,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(0,1) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(2,0);
        const T ctheta1 = this->data(1,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = this->data(2,1) / ctheta2;
    const T ctheta1 = this->data(1,1) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(0,2) / ctheta2;
    const T ctheta3 = this->data(0,0) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body213() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(-this->data(1,2));

    if (std::abs(this->data(1,2) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(0,1);
        const T ctheta1 = -this->data(2,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(1,2) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(0,1);
        const T ctheta1 = this->data(2,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = this->data(0,2) / ctheta2;
    const T ctheta1 = this->data(2,2) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(1,0) / ctheta2;
    const T ctheta3 = this->data(1,1) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body321() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(-this->data(2,0));

    if (std::abs(this->data(2,0) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(1,2);
        const T ctheta1 = -this->data(0,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(2,0) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(1,2);
        const T ctheta1 = this->data(0,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = this->data(1,0) / ctheta2;
    const T ctheta1 = this->data(0,0) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(2,1) / ctheta2;
    const T ctheta3 = this->data(2,2) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body121() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(0,0));

    if (std::abs(this->data(0,0) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(1,2);
        const T ctheta1 = this->data(2,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(0,0) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(1,2);
        const T ctheta1 = -this->data(2,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(1,0) / stheta2;
    const T ctheta1 = -this->data(2,0) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(0,1) / stheta2;
    const T ctheta3 = this->data(0,2) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body131() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(0,0));

    if (std::abs(this->data(0,0) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(2,1);
        const T ctheta1 = this->data(1,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(0,0) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(2,1);
        const T ctheta1 = -this->data(1,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(2,0) / stheta2;
    const T ctheta1 = this->data(1,0) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(0,2) / stheta2;
    const T ctheta3 = -this->data(0,1) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body212() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(1,1));

    if (std::abs(this->data(1,1) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(0,2);
        const T ctheta1 = this->data(2,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(1,1) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(0,2);
        const T ctheta1 = -this->data(2,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(0,1) / stheta2;
    const T ctheta1 = this->data(2,1) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(1,0) / stheta2;
    const T ctheta3 = -this->data(1,2) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body232() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(1,1));

    if (std::abs(this->data(1,1) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(2,0);
        const T ctheta1 = this->data(0,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(1,1) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(2,0);
        const T ctheta1 = -this->data(0,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(2,1) / stheta2;
    const T ctheta1 = -this->data(0,1) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(1,2) / stheta2;
    const T ctheta3 = this->data(1,0) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body313() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(2,2));

    if (std::abs(this->data(2,2) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(0,1);
        const T ctheta1 = this->data(1,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(2,2) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(0,1);
        const T ctheta1 = -this->data(1,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(0,2) / stheta2;
    const T ctheta1 = -this->data(1,2) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(2,0) / stheta2;
    const T ctheta3 = this->data(2,1) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_body323() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(2,2));

    if (std::abs(this->data(2,2) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(1,0);
        const T ctheta1 = this->data(0,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(2,2) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(1,0);
        const T ctheta1 = -this->data(0,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(1,2) / stheta2;
    const T ctheta1 = this->data(0,2) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(2,1) / stheta2;
    const T ctheta3 = -this->data(2,0) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space123() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(-this->data(2,0));

    if (std::abs(this->data(2,0) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(0,1);
        const T ctheta1 = -this->data(0,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(2,0) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(0,1);
        const T ctheta1 = this->data(0,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = this->data(2,1) / ctheta2;
    const T ctheta1 = this->data(2,2) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(1,0) / ctheta2;
    const T ctheta3 = this->data(0,0) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space231() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(-this->data(0,1));

    if (std::abs(this->data(0,1) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(1,2);
        const T ctheta1 = -this->data(1,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(0,1) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(1,2);
        const T ctheta1 = this->data(1,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = this->data(0,2) / ctheta2;
    const T ctheta1 = this->data(0,0) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(2,1) / ctheta2;
    const T ctheta3 = this->data(1,1) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space312() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(-this->data(1,2));

    if (std::abs(this->data(1,2) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(2,0);
        const T ctheta1 = -this->data(2,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(1,2) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(2,0);
        const T ctheta1 = this->data(2,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = this->data(1,0) / ctheta2;
    const T ctheta1 = this->data(1,1) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(0,2) / ctheta2;
    const T ctheta3 = this->data(2,2) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space132() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(this->data(1,0));

    if (std::abs(this->data(1,0) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(0,2);
        const T ctheta1 = -this->data(0,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(1,0) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(0,2);
        const T ctheta1 = this->data(0,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = -this->data(1,2) / ctheta2;
    const T ctheta1 = this->data(1,1) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = -this->data(2,0) / ctheta2;
    const T ctheta3 = this->data(0,0) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


template <typename T>
std::array<T, 3> SO::to_eulerangles_space213() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(this->data(2,1));

    if (std::abs(this->data(2,1) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(1,0);
        const T ctheta1 = -this->data(1,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(2,1) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(1,0);
        const T ctheta1 = this->data(1,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = -this->data(2,0) / ctheta2;
    const T ctheta1 = this->data(2,2) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = -this->data(0,1) / ctheta2;
    const T ctheta3 = this->data(1,1) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space321() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::asin(this->data(0,2));

    if (std::abs(this->data(0,2) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(2,1);
        const T ctheta1 = -this->data(2,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(0,2) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(2,1);
        const T ctheta1 = this->data(2,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = std::cos(theta2);

    const T stheta1 = -this->data(0,1) / ctheta2;
    const T ctheta1 = this->data(0,0) / ctheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = -this->data(1,2) / ctheta2;
    const T ctheta3 = this->data(2,2) / ctheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space121() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(0,0));

    if (std::abs(this->data(0,0) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(2,1);
        const T ctheta1 = this->data(2,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(0,0) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(2,1);
        const T ctheta1 = -this->data(2,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(0,1) / stheta2;
    const T ctheta1 = this->data(0,2) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(1,0) / stheta2;
    const T ctheta3 = -this->data(2,0) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space131() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(0,0));

    if (std::abs(this->data(0,0) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(1,2);
        const T ctheta1 = this->data(1,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(0,0) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(1,2);
        const T ctheta1 = -this->data(1,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(0,2) / stheta2;
    const T ctheta1 = -this->data(0,1) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(2,0) / stheta2;
    const T ctheta3 = this->data(1,0) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space212() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(1,1));

    if (std::abs(this->data(1,1) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(2,0);
        const T ctheta1 = this->data(2,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(1,1) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(2,0);
        const T ctheta1 = -this->data(2,2);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(1,0) / stheta2;
    const T ctheta1 = -this->data(1,2) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(0,1) / stheta2;
    const T ctheta3 = this->data(2,1) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space232() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(1,1));

    if (std::abs(this->data(1,1) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(0,2);
        const T ctheta1 = this->data(0,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(1,1) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(0,2);
        const T ctheta1 = -this->data(0,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(1,2) / stheta2;
    const T ctheta1 = this->data(1,0) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(2,1) / stheta2;
    const T ctheta3 = -this->data(0,1) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space313() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(2,2));

    if (std::abs(this->data(2,2) - 1.0) < 1e-14)
    {
        const T stheta1 = this->data(1,0);
        const T ctheta1 = this->data(1,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(2,2) + 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(1,0);
        const T ctheta1 = -this->data(1,1);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(2,0) / stheta2;
    const T ctheta1 = this->data(2,1) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(0,2) / stheta2;
    const T ctheta3 = -this->data(1,2) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> SO::to_eulerangles_space323() const
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] this->data The direction co-sine matrix, \f$this->data\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    if (this->_shape != 3)
    {
        throw std::domain_error("SO::to_eulerangles: Expected input shape 3. Got " + std::to_string(this->_shape) + ".");
    }

    const T theta2 = std::acos(this->data(2,2));

    if (std::abs(this->data(2,2) - 1.0) < 1e-14)
    {
        const T stheta1 = -this->data(0,1);
        const T ctheta1 = this->data(0,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    if (std::abs(this->data(2,2) + 1.0) < 1e-14)
    {
        const T stheta1 = this->data(0,1);
        const T ctheta1 = -this->data(0,0);
        const T theta1 = std::atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = std::sin(theta2);

    const T stheta1 = this->data(2,1) / stheta2;
    const T ctheta1 = -this->data(2,0) / stheta2;
    const T theta1 = std::atan2(stheta1, ctheta1);

    const T stheta3 = this->data(1,2) / stheta2;
    const T ctheta3 = this->data(0,2) / stheta2;
    const T theta3 = std::atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

}

#endif
