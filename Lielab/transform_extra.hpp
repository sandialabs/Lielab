#ifndef _LIELAB_TRANSFORM_HPP
#define _LIELAB_TRANSFORM_HPP

// Standard C++ includes
#include <array>
#include <cmath>
#include <exception>
#include <tuple>

// lielab includes
#include "domain.hpp"

// 3rd party includes
#include <Eigen/Core>

using namespace std; // asin, abs, atan2, cos

namespace Lielab
{
/*!
* The transform namespace. Houses all functions and classes related to coordinate transformations.
*/
namespace transform
{
Lielab::domain::SU dcm_to_quaternion(const Lielab::domain::SO &dcm)
{
    /*!
        * Transforms a Quaternion to a Direction Cosine object.
        *
        * @param[in] dcm A direction cosine as an SO object of shape 3.
        * @param[out] q A quaternion as an SU object of shape 2.
        */

    if (dcm.shape != 3)
    {
        throw std::domain_error("dcm_to_quaternion: Expected input shape 3. Got " + std::to_string(dcm.shape) + ".");
    }

    const double e[] = {0.5*std::sqrt(1.0 + dcm(0,0) + dcm(1,1) + dcm(2,2)),
                        0.5*std::sqrt(1.0 + dcm(0,0) - dcm(1,1) - dcm(2,2)),
                        0.5*std::sqrt(1.0 - dcm(0,0) + dcm(1,1) - dcm(2,2)),
                        0.5*std::sqrt(1.0 - dcm(0,0) - dcm(1,1) + dcm(2,2))};
    
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
        return Lielab::domain::SU::from_quaternion(e[0],
                                                (dcm(2,1) - dcm(1,2))/(4.0*e[0]),
                                                (dcm(0,2) - dcm(2,0))/(4.0*e[0]),
                                                (dcm(1,0) - dcm(0,1))/(4.0*e[0]));
    }
    else if (ind == 1)
    {
        return Lielab::domain::SU::from_quaternion((dcm(2,1) - dcm(1,2))/(4.0*e[1]),
                                                e[1],
                                                (dcm(1,0) + dcm(0,1))/(4.0*e[1]),
                                                (dcm(0,2) + dcm(2,0))/(4.0*e[1]));
    }
    else if (ind == 2)
    {
        return Lielab::domain::SU::from_quaternion((dcm(0,2) - dcm(2,0))/(4.0*e[2]),
                                                (dcm(1,0) + dcm(0,1))/(4.0*e[2]),
                                                e[2],
                                                (dcm(2,1) + dcm(1,2))/(4.0*e[2]));
    }
    else
    {
        // ind == 3 should be the only other option by this point.
        assert(ind == 3);

        return Lielab::domain::SU::from_quaternion((dcm(1,0) - dcm(0,1))/(4.0*e[3]),
                                                (dcm(0,2) + dcm(2,0))/(4.0*e[3]),
                                                (dcm(2,1) + dcm(1,2))/(4.0*e[3]),
                                                e[3]);
    }
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody123(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = std::asin(dcm(0,2));

    if ((abs(dcm(0,2) - 1.0) < 1e-14) || abs(dcm(0,2) + 1.0) < 1e-14)
    {
        const T stheta1 = dcm(1,0);
        const T ctheta1 = dcm(1,1);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = -dcm(1,2) / ctheta2;
    const T ctheta1 = dcm(2,2) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = -dcm(0,1) / ctheta2;
    const T ctheta3 = dcm(0,0) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody231(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(dcm(1,0));

    if ((abs(dcm(1,0) - 1) < 1e-14) || abs(dcm(1,0) + 1) < 1e-14)
    {
        const T stheta1 = dcm(0,2);
        const T ctheta1 = -dcm(0,1);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = -dcm(2,0) / ctheta2;
    const T ctheta1 = dcm(0,0) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = -dcm(1,2) / ctheta2;
    const T ctheta3 = dcm(1,1) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody312(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(dcm(2,1));

    if ((abs(dcm(2,1) - 1) < 1e-14) || abs(dcm(2,1) + 1) < 1e-14)
    {
        const T stheta1 = dcm(0,2);
        const T ctheta1 = dcm(0,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = -dcm(0,1) / ctheta2;
    const T ctheta1 = dcm(1,1) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = -dcm(2,0) / ctheta2;
    const T ctheta3 = dcm(2,2) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody132(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(-dcm(0,1));

    if ((std::abs(dcm(0,1) - 1) < 1e-14) || abs(dcm(0,1) + 1) < 1e-14)
    {
        const T stheta1 = -dcm(1,2);
        const T ctheta1 = dcm(1,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = dcm(2,1) / ctheta2;
    const T ctheta1 = dcm(1,1) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(0,2) / ctheta2;
    const T ctheta3 = dcm(0,0) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody213(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(-dcm(1,2));

    if ((abs(dcm(1,2) - 1) < 1e-14) || abs(dcm(1,2) + 1) < 1e-14)
    {
        const T stheta1 = -dcm(0,1);
        const T ctheta1 = dcm(0,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = dcm(0,2) / ctheta2;
    const T ctheta1 = dcm(2,2) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(1,0) / ctheta2;
    const T ctheta3 = dcm(1,1) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody321(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(-dcm(2,0));

    if ((abs(dcm(2,0) - 1) < 1e-14) || abs(dcm(2,0) + 1) < 1e-14)
    {
        const T stheta1 = -dcm(0,1);
        const T ctheta1 = dcm(0,2);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = dcm(1,0) / ctheta2;
    const T ctheta1 = dcm(0,0) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(2,1) / ctheta2;
    const T ctheta3 = dcm(2,2) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody121(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(0,0));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = -dcm(1,2);
        const T ctheta1 = dcm(1,1);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(1,0) / stheta2;
    const T ctheta1 = -dcm(2,0) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(0,1) / stheta2;
    const T ctheta3 = dcm(0,2) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody131(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(0,0));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = -dcm(1,2);
        const T ctheta1 = dcm(1,1);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(2,0) / stheta2;
    const T ctheta1 = dcm(1,0) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(0,2) / stheta2;
    const T ctheta3 = -dcm(0,1) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody212(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(1,1));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = dcm(0,2);
        const T ctheta1 = dcm(0,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(0,1) / stheta2;
    const T ctheta1 = dcm(2,1) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(1,0) / stheta2;
    const T ctheta3 = -dcm(1,2) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody232(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(1,1));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = dcm(0,2);
        const T ctheta1 = dcm(0,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(2,1) / stheta2;
    const T ctheta1 = -dcm(0,1) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(1,2) / stheta2;
    const T ctheta3 = dcm(1,0) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody313(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(2,2));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = -dcm(0,1);
        const T ctheta1 = dcm(0,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(0,2) / stheta2;
    const T ctheta1 = -dcm(1,2) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(2,0) / stheta2;
    const T ctheta3 = dcm(2,1) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglebody323(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(2,2));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = -dcm(0,1);
        const T ctheta1 = dcm(0,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(1,2) / stheta2;
    const T ctheta1 = dcm(0,2) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(2,1) / stheta2;
    const T ctheta3 = -dcm(2,0) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace123(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(-dcm(2,0));

    if ((abs(dcm(2,0) - 1.0) < 1e-14) || abs(dcm(2,0) + 1.0) < 1e-14)
    {
        const T stheta1 = dcm(0,1);
        const T ctheta1 = dcm(0,2);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = dcm(2,1) / ctheta2;
    const T ctheta1 = dcm(2,2) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(1,0) / ctheta2;
    const T ctheta3 = dcm(0,0) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace231(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(-dcm(0,1));

    if ((abs(dcm(0,1) - 1) < 1e-14) || abs(dcm(0,1) + 1) < 1e-14)
    {
        const T stheta1 = dcm(1,2);
        const T ctheta1 = dcm(1,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = dcm(0,2) / ctheta2;
    const T ctheta1 = dcm(0,0) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(2,1) / ctheta2;
    const T ctheta3 = dcm(1,1) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace312(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(-dcm(1,2));

    if ((abs(dcm(1,2) - 1) < 1e-14) || abs(dcm(1,2) + 1) < 1e-14)
    {
        const T stheta1 = dcm(2,0);
        const T ctheta1 = dcm(2,1);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = dcm(1,0) / ctheta2;
    const T ctheta1 = dcm(1,1) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(0,2) / ctheta2;
    const T ctheta3 = dcm(2,2) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace132(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(dcm(1,0));

    if ((abs(dcm(1,0) - 1) < 1e-14) || abs(dcm(1,0) + 1) < 1e-14)
    {
        const T stheta1 = dcm(0,2);
        const T ctheta1 = -dcm(0,1);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = -dcm(1,2) / ctheta2;
    const T ctheta1 = dcm(1,1) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = -dcm(2,0) / ctheta2;
    const T ctheta3 = dcm(0,0) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}


template <typename T>
std::array<T, 3> dcm_to_eanglespace213(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(dcm(2,1));

    if ((abs(dcm(2,1) - 1) < 1e-14) || abs(dcm(2,1) + 1) < 1e-14)
    {
        const T stheta1 = dcm(1,0);
        const T ctheta1 = -dcm(1,2);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = -dcm(2,0) / ctheta2;
    const T ctheta1 = dcm(2,2) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = -dcm(0,1) / ctheta2;
    const T ctheta3 = dcm(1,1) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace321(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = asin(dcm(0,2));

    if ((abs(dcm(0,2) - 1) < 1e-14) || abs(dcm(0,2) + 1) < 1e-14)
    {
        const T stheta1 = dcm(2,1);
        const T ctheta1 = -dcm(2,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T ctheta2 = cos(theta2);

    const T stheta1 = -dcm(0,1) / ctheta2;
    const T ctheta1 = dcm(0,0) / ctheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = -dcm(1,2) / ctheta2;
    const T ctheta3 = dcm(2,2) / ctheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace121(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(0,0));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = -dcm(1,2);
        const T ctheta1 = dcm(1,1);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(0,1) / stheta2;
    const T ctheta1 = dcm(0,2) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(1,0) / stheta2;
    const T ctheta3 = -dcm(2,0) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace131(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(0,0));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = -dcm(1,2);
        const T ctheta1 = dcm(1,1);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(0,2) / stheta2;
    const T ctheta1 = -dcm(0,1) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(2,0) / stheta2;
    const T ctheta3 = dcm(1,0) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace212(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(1,1));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = dcm(0,2);
        const T ctheta1 = dcm(0,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(1,0) / stheta2;
    const T ctheta1 = -dcm(1,2) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(0,1) / stheta2;
    const T ctheta3 = dcm(2,1) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace232(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(1,1));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = dcm(0,2);
        const T ctheta1 = dcm(0,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(1,2) / stheta2;
    const T ctheta1 = dcm(1,0) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(2,1) / stheta2;
    const T ctheta3 = -dcm(0,1) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace313(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(2,2));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = -dcm(0,1);
        const T ctheta1 = dcm(0,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(2,0) / stheta2;
    const T ctheta1 = dcm(2,1) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(0,2) / stheta2;
    const T ctheta3 = -dcm(1,2) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

template <typename T>
std::array<T, 3> dcm_to_eanglespace323(const Lielab::domain::SO &dcm)
{
    /*! \f{equation*}{ SO(3) \rightarrow (\mathbb{R}, \mathbb{R}, \mathbb{R}) \f}
    * 
    * @param[in] dcm The direction co-sine matrix, \f$dcm\f$.
    * @param[out] theta1 Rotation "1" in radians, \f$\theta_1\f$.
    * @param[out] theta2 Rotation "2" in radians, \f$\theta_2\f$.
    * @param[out] theta3 Rotation "3" in radians, \f$\theta_3\f$.
    */

    const T theta2 = acos(dcm(2,2));

    if (abs(theta2) < 1e-14)
    {
        const T stheta1 = -dcm(0,1);
        const T ctheta1 = dcm(0,0);
        const T theta1 = atan2(stheta1, ctheta1);
        return {theta1, theta2, 0.0};
    }

    const T stheta2 = sin(theta2);

    const T stheta1 = dcm(2,1) / stheta2;
    const T ctheta1 = -dcm(2,0) / stheta2;
    const T theta1 = atan2(stheta1, ctheta1);

    const T stheta3 = dcm(1,2) / stheta2;
    const T ctheta3 = dcm(0,2) / stheta2;
    const T theta3 = atan2(stheta3, ctheta3);

    return {theta1, theta2, theta3};
}

Eigen::VectorXd dcm_to_gibbs(const Lielab::domain::SO &dcm)
{
    /*!
        * Transforms a direction co-sine matrics into a Gibbs (Rodrigues) vector.
        * @param[in] dcm The direction co-sine matrix.
        * @param[out] gibbs The rotation in Gibbs representation.
        */

    Eigen::VectorXd gibbs(3);
    const double mul = 1.0/(1.0 + dcm(0,0) + dcm(1,1) + dcm(2,2));
    
    gibbs(0) = dcm(1,2) - dcm(2,1);
    gibbs(1) = dcm(2,0) - dcm(0,2);
    gibbs(2) = dcm(0,1) - dcm(1,0);
    gibbs = gibbs*mul;

    return gibbs;
}

Lielab::domain::SO quaternion_to_dcm(const Lielab::domain::SU &q)
{
    /*!
        * Transforms a Quaternion to a Direction Cosine object.
        *
        * @param[in] q A quaternion as an SU object of shape 2.
        * @param[out] dcm A direction cosine as an SO object of shape 3.
        */

    if (q.shape != 2)
    {
        throw std::domain_error("quaternion_to_dcm: Expected input shape 2. Got " + std::to_string(q.shape) + ".");
    }

    // TODO: Fix for new SU serialization

    // const Eigen::VectorXd q_real = q.serialize();
    // const double e0 = q_real(0);
    // const double e1 = q_real(1);
    // const double e2 = q_real(2);
    // const double e3 = q_real(3);

    const double e0 = (std::real(q._data(0,0)) + std::real(q._data(1,1))) / 2.0;
    const double e1 = (std::imag(q._data(0,0)) - std::imag(q._data(1,1))) / 2.0;
    const double e2 = (std::real(q._data(1,0)) - std::real(q._data(0,1))) / 2.0;
    const double e3 = (std::imag(q._data(1,0)) + std::imag(q._data(0,1))) / 2.0;

    Lielab::domain::SO dcm(3);
    dcm(0,0) = 1.0 - 2.0*e2*e2 - 2.0*e3*e3;
    dcm(0,1) = 2.0*e1*e2 - 2.0*e3*e0;
    dcm(0,2) = 2.0*e1*e3 + 2.0*e2*e0;
    dcm(1,0) = 2.0*e1*e2 + 2.0*e3*e0;
    dcm(1,1) = 1.0 - 2.0*e1*e1 - 2.0*e3*e3;
    dcm(1,2) = 2.0*e2*e3 - 2.0*e1*e0;
    dcm(2,0) = 2.0*e1*e3 - 2.0*e2*e0;
    dcm(2,1) = 2.0*e2*e3 + 2.0*e1*e0;
    dcm(2,2) = 1.0 - 2.0*e1*e1 - 2.0*e2*e2;

    return dcm;
}
}
}

#endif