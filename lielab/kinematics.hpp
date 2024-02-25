#ifndef _LIELAB_KINEMATICS_H
#define _LIELAB_KINEMATICS_H

#include <Eigen/Core>
#include <cmath>

namespace lielab
{
    /*!
    * The kinematics namespace. Houses all functions and classes related to kinematics and their differential equations.
    */
    namespace kinematics
    {
        Eigen::MatrixXd dcm(Eigen::VectorXd angular_rates, Eigen::MatrixXd dcm)
        {
            /*!
            * Computes the rate-of-change of a direction co-sine matrix.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] dcm The current direction co-sine value.
            * @param[out] ddcm The rate-of-change of the direction co-sine matrix.
            */
            Eigen::MatrixXd ddcm(3,3);
            Eigen::VectorXd w = angular_rates;

            ddcm(0,0) = w(2)*dcm(1,0) - w(1)*dcm(2,0);
            ddcm(0,1) = w(2)*dcm(1,1) - w(1)*dcm(2,1);
            ddcm(0,2) = w(2)*dcm(1,2) - w(1)*dcm(2,2);
            ddcm(1,0) = w(0)*dcm(2,0) - w(2)*dcm(0,0);
            ddcm(1,1) = w(0)*dcm(2,1) - w(2)*dcm(0,1);
            ddcm(1,2) = w(0)*dcm(2,2) - w(2)*dcm(0,2);
            ddcm(2,0) = w(1)*dcm(0,0) - w(0)*dcm(1,0);
            ddcm(2,1) = w(1)*dcm(0,1) - w(0)*dcm(1,1);
            ddcm(2,2) = w(1)*dcm(0,2) - w(0)*dcm(1,2);

            return ddcm;
        }

        Eigen::VectorXd quaternions(Eigen::VectorXd angular_rates, lielab::domain::SU q)
        {
            /*!
            * Computes the rate-of-change of quaternions.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] q The current quaternion value.
            * @param[out] dq The rate-of-change of the quaternion.
            */

            Eigen::VectorXd qvec = q.serialize();
            Eigen::VectorXd dq(4);
            Eigen::VectorXd w = angular_rates;

            dq(0) = 1.0/2.0*(w(0)*qvec(3) - w(1)*qvec(2) + w(2)*qvec(1));
            dq(1) = 1.0/2.0*(w(0)*qvec(2) + w(1)*qvec(3) - w(2)*qvec(0));
            dq(2) = 1.0/2.0*(-w(0)*qvec(1) + w(1)*qvec(0) + w(2)*qvec(3));
            dq(3) = 1.0/2.0*(-w(0)*qvec(0) - w(1)*qvec(1) - w(2)*qvec(2));

            // TODO: Literature seems to be inconsistent with quaternion definitions
            //       and SU. Sometimes its transposed in the literature but
            //       I don't believe thats correct.
            // dq(0) = 1.0/2.0*(-w(0)*q(1) - w(1)*q(2) - w(2)*q(3));
            // dq(1) = 1.0/2.0*(w(0)*q(0) + w(2)*q(2) - w(1)*q(3));
            // dq(2) = 1.0/2.0*(w(1)*q(0) - w(2)*q(1) + w(0)*q(3));
            // dq(3) = 1.0/2.0*(w(2)*q(0) + w(1)*q(1) - w(0)*q(2));

            return dq;
        }

        Eigen::VectorXd eanglespace123(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 123 rotation.
            * @param[in] angular_rates The angular rate in radians per second
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = w(0) + (w(2)*std::cos(angle1)*std::sin(angle2))/std::cos(angle2) + (w(1)*std::sin(angle1)*std::sin(angle2))/std::cos(angle2);
            out(1) = w(1)*std::cos(angle1) - w(2)*std::sin(angle1);
            out(2) = (w(2)*std::cos(angle1))/std::cos(angle2) + (w(1)*std::sin(angle1))/std::cos(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace231(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 231 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(2)*std::sin(angle1) + w(0)*std::cos(angle1))*(std::sin(angle2)/std::cos(angle2)) + w(1);
            out(1) = w(2)*std::cos(angle1) - w(0)*std::sin(angle1);
            out(2) = (w(2)*std::sin(angle1) + w(0)*std::cos(angle1))/std::cos(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace312(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 312 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(0)*std::sin(angle1) + w(1)*std::cos(angle1))*(std::sin(angle2)/std::cos(angle2)) + w(2);
            out(1) = w(0)*std::cos(angle1) - w(1)*std::sin(angle1);
            out(2) = (w(0)*std::sin(angle1) + w(1)*std::cos(angle1))/std::cos(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace132(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 132 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(2)*std::sin(angle1) - w(1)*std::cos(angle1))*(std::sin(angle2)/std::cos(angle2)) + w(0);
            out(1) = w(2)*std::cos(angle1) + w(1)*std::sin(angle1);
            out(2) = (-w(2)*std::sin(angle1) + w(1)*std::cos(angle1))/std::cos(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace213(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 213 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(0)*std::sin(angle1) - w(2)*std::cos(angle1))*(std::sin(angle2)/std::cos(angle2)) + w(2);
            out(1) = w(0)*std::cos(angle1) - w(2)*std::sin(angle1);
            out(2) = (-w(0)*std::sin(angle1) + w(2)*std::cos(angle1))/std::cos(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace321(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 321 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(1)*std::sin(angle1) - w(0)*std::cos(angle1))*(std::sin(angle2)/std::cos(angle2)) + w(2);
            out(1) = w(1)*std::cos(angle1) + w(0)*std::sin(angle1);
            out(2) = (-w(1)*std::sin(angle1) + w(0)*std::cos(angle1))/std::cos(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace121(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 121 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = -(w(1)*std::sin(angle1) + w(2)*std::cos(angle1))*(std::cos(angle2)/std::sin(angle2)) + w(0);
            out(1) = w(1)*std::cos(angle1) - w(2)*std::sin(angle1);
            out(2) = (w(1)*std::sin(angle1) + w(2)*std::cos(angle1))/std::sin(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace131(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 131 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (-w(2)*std::sin(angle1) + w(1)*std::cos(angle1))*(std::cos(angle2)/std::sin(angle2)) + w(0);
            out(1) = w(2)*std::cos(angle1) + w(1)*std::sin(angle1);
            out(2) = (w(2)*std::sin(angle1) - w(1)*std::cos(angle1))/std::sin(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace212(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 212 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (-w(0)*std::sin(angle1) + w(2)*std::cos(angle1))*(std::cos(angle2)/std::sin(angle2)) + w(1);
            out(1) = w(0)*std::cos(angle1) + w(2)*std::sin(angle1);
            out(2) = (w(0)*std::sin(angle1) - w(2)*std::cos(angle1))/std::sin(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace232(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 232 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = -(w(2)*std::sin(angle1) + w(0)*std::cos(angle1))*(std::cos(angle2)/std::sin(angle2)) + w(1);
            out(1) = w(2)*std::cos(angle1) + w(0)*std::sin(angle1);
            out(2) = (w(2)*std::sin(angle1) + w(0)*std::cos(angle1))/std::sin(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace313(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 313 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = -(w(0)*std::sin(angle1) + w(1)*std::cos(angle1))*(std::cos(angle2)/std::sin(angle2)) + w(2);
            out(1) = w(0)*std::cos(angle1) + w(1)*std::sin(angle1);
            out(2) = (w(0)*std::sin(angle1) + w(1)*std::cos(angle1))/std::sin(angle2);

            return out;
        }

        Eigen::VectorXd eanglespace323(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler space 323 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the euler angle.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (-w(1)*std::sin(angle1) + w(0)*std::cos(angle1))*(std::cos(angle2)/std::sin(angle2)) + w(2);
            out(1) = w(1)*std::cos(angle1) + w(0)*std::sin(angle1);
            out(2) = (w(1)*std::sin(angle1) - w(0)*std::cos(angle1))/std::sin(angle2);

            return out;
        }

        Eigen::VectorXd eanglebody123(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 123 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(0)*std::cos(angle3) - w(1)*std::sin(angle3)) / std::cos(angle2);
            out(1) = w(0)*std::sin(angle3) + w(1)*std::cos(angle3);
            out(2) = (w(1)*std::sin(angle3) - w(0)*std::cos(angle3))*(std::sin(angle2) / std::cos(angle2)) + w(2);

            return out;
        }

        Eigen::VectorXd eanglebody231(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 231 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(1)*std::cos(angle3) - w(2)*std::sin(angle3)) / std::cos(angle2);
            out(1) = w(1)*std::sin(angle3) + w(2)*std::cos(angle3);
            out(2) = w(0) + (w(2)*std::sin(angle3) - w(1)*std::cos(angle3))*(std::sin(angle2)/std::cos(angle2));

            return out;
        }

        Eigen::VectorXd eanglebody312(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 312 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(2)*std::cos(angle3) - w(0)*std::sin(angle3)) / std::cos(angle2);
            out(1) = w(0)*std::cos(angle3) + w(2)*std::sin(angle3);
            out(2) = w(1) + (w(0)*std::sin(angle3) - w(2)*std::cos(angle3))*(std::sin(angle2)/std::cos(angle2));

            return out;
        }

        Eigen::VectorXd eanglebody132(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 132 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(0)*std::cos(angle3) + w(2)*std::sin(angle3)) / std::cos(angle2);
            out(1) = -w(0)*std::sin(angle3) + w(2)*std::cos(angle3);
            out(2) = w(1) + (w(2)*std::sin(angle3) + w(0)*std::cos(angle3))*(std::sin(angle2)/std::cos(angle2));

            return out;
        }

        Eigen::VectorXd eanglebody213(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 213 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(1)*std::cos(angle3) + w(0)*std::sin(angle3)) / std::cos(angle2);
            out(1) = -w(1)*std::sin(angle3) + w(0)*std::cos(angle3);
            out(2) = w(2) + (w(0)*std::sin(angle3) + w(1)*std::cos(angle3))*(std::sin(angle2)/std::cos(angle2));

            return out;
        }

        Eigen::VectorXd eanglebody321(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 321 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(2)*std::cos(angle3) + w(1)*std::sin(angle3)) / std::cos(angle2);
            out(1) = -w(2)*std::sin(angle3) + w(1)*std::cos(angle3);
            out(2) = w(2) + (w(1)*std::sin(angle3) + w(2)*std::cos(angle3))*(std::sin(angle2)/std::cos(angle2));

            return out;
        }

        Eigen::VectorXd eanglebody121(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 121 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(2)*std::cos(angle3) + w(1)*std::sin(angle3)) / std::sin(angle2);
            out(1) = -w(2)*std::sin(angle3) + w(1)*std::cos(angle3);
            out(2) = w(2) - (w(1)*std::sin(angle3) + w(2)*std::cos(angle3))*(std::cos(angle2)/std::sin(angle2));

            return out;
        }

        Eigen::VectorXd eanglebody131(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 131 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (-w(1)*std::cos(angle3) + w(2)*std::sin(angle3)) / std::sin(angle2);
            out(1) = w(1)*std::sin(angle3) + w(2)*std::cos(angle3);
            out(2) = w(0) + (-w(2)*std::sin(angle3) + w(1)*std::cos(angle3))*(std::cos(angle2)/std::sin(angle2));

            return out;
        }

        Eigen::VectorXd eanglebody212(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 212 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(2)*std::cos(angle3) + w(0)*std::sin(angle3)) / std::sin(angle2);
            out(1) = w(2)*std::sin(angle3) + w(0)*std::cos(angle3);
            out(2) = w(1) + (-w(0)*std::sin(angle3) + w(2)*std::cos(angle3))*(std::cos(angle2)/std::sin(angle2));

            return out;
        }

        Eigen::VectorXd eanglebody232(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 232 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(0)*std::cos(angle3) + w(2)*std::sin(angle3)) / std::sin(angle2);
            out(1) = -w(0)*std::sin(angle3) + w(2)*std::cos(angle3);
            out(2) = w(1) - (w(2)*std::sin(angle3) + w(0)*std::cos(angle3))*(std::cos(angle2)/std::sin(angle2));

            return out;
        }

        Eigen::VectorXd eanglebody313(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 313 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(1)*std::cos(angle3) + w(0)*std::sin(angle3)) / std::sin(angle2);
            out(1) = -w(1)*std::sin(angle3) + w(0)*std::cos(angle3);
            out(2) = w(2) - (w(0)*std::sin(angle3) + w(1)*std::cos(angle3))*(std::cos(angle2)/std::sin(angle2));

            return out;
        }

        Eigen::VectorXd eanglebody323(Eigen::VectorXd angular_rates, double angle1, double angle2, double angle3)
        {
            /*!
            * Computes the rate-of-change of an euler body 323 rotation.
            * @param[in] angular_rates The angular rate in radians per second.
            * @param[in] angle1 First rotation in radians.
            * @param[in] angle2 Second rotation in radians.
            * @param[in] angle3 Third rotation in radians.
            * @param[out] out The rate-of-change of the rotation.
            */
            Eigen::VectorXd out(3);
            Eigen::VectorXd w = angular_rates;

            out(0) = (w(0)*std::cos(angle3) + w(1)*std::sin(angle3)) / std::sin(angle2);
            out(1) = w(0)*std::sin(angle3) + w(1)*std::cos(angle3);
            out(2) = w(2) + (-w(1)*std::sin(angle3) + w(0)*std::cos(angle3))*(std::cos(angle2)/std::sin(angle2));

            return out;
        }
    }
}

#endif
