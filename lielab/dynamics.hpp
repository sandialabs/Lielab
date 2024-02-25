#ifndef _LIELAB_DYNAMICS_H
#define _LIELAB_DYNAMICS_H

#include <Eigen/Core>
#include <cmath>
#include <variant>

using std::sin;
using std::cos;
using std::pow;

#include "domain"
#include "topos.hpp"

namespace lielab
{
    /*!
    * The dynamics namespace. Houses all functions and classes related to dynamics and their differential equations.
    */
    namespace dynamics
    {
        typedef std::function<lielab::domain::halie(double, lielab::domain::hmlie)> vectorfield;

        vectorfield vfex1()
        {
            /*!
            * Constructs and returns eom vfex1.
            */

            vectorfield out = [](const double t, const lielab::domain::hmlie & M)
            {
                lielab::domain::SO y = std::get<lielab::domain::SO>(M.space[0]);
                lielab::domain::so dy{std::pow(t, 2), 1.0, -t};
                return lielab::domain::halie{dy};
            };

            return out;
        }

        vectorfield vfex2()
        {
            vectorfield out = [](const double t, const lielab::domain::hmlie & M)
            {
                lielab::domain::RN y = std::get<lielab::domain::RN>(M.space[0]);

                const double sigma = 10.0;
                const double rho = 28.0;
                const double b1 = 8.0;
                const double b2 = 3.0;
                const double beta = b1/b2;

                lielab::domain::rn dx(4);
                dx(0) = -beta*y(0) + y(1)*y(2);
                dx(1) = -sigma*y(1) + sigma*y(2);
                dx(2) = -y(0)*y(1) + rho*y(1) - y(2);

                return lielab::domain::halie{dx};
            };

            return out;
        }

        // vectorfield vfex3()
        // {
        //     vectorfield out = [](const double t, const lielab::domain::hmlie & M)
        //     {

        //     };

        //     return out;
        // }

        vectorfield rigidbody_real_dcm(Eigen::Matrix3d inertia)
        {
            Eigen::Matrix3d invinertia = inertia.inverse();

            vectorfield out = [invinertia, inertia](const double _t, const lielab::domain::hmlie & _y)
            {
                lielab::domain::RN inp = std::get<lielab::domain::RN>(_y.space[0]);
                lielab::domain::rn dx(13);
                Eigen::Vector3d omega;

                double Ixx = inertia(0,0);
                double Iyy = inertia(1,1);
                double Izz = inertia(2,2);

                double w1 = inp(0);
                double w2 = inp(1);
                double w3 = inp(2);
                double c11 = inp(3);
                double c12 = inp(4);
                double c13 = inp(5);
                double c21 = inp(6);
                double c22 = inp(7);
                double c23 = inp(8);
                double c31 = inp(9);
                double c32 = inp(10);
                double c33 = inp(11);

                omega(0) = w1;
                omega(1) = w2;
                omega(2) = w3;

                Eigen::Vector3d pdot;
                // pdot = invinertia*(-omega.cross(omega.transpose()*inertia));

                double f = 0.0;
                if (_t >= 8.0 && _t <= 10.0)
                {
                    f = 1.0;
                }

                pdot(0) = -omega(1)*omega(2)*(Izz - Iyy)/Ixx + f/Ixx;
                pdot(1) = -omega(0)*omega(2)*(Ixx - Izz)/Iyy;
                pdot(2) = -omega(1)*omega(0)*(Iyy - Ixx)/Izz;

                dx(0) = pdot(0);
                dx(1) = pdot(1);
                dx(2) = pdot(2);
                dx(3) = w3*c21 - w2*c31;
                dx(4) = w3*c22 - w2*c32;
                dx(5) = w3*c23 - w2*c33;
                dx(6) = w1*c31 - w3*c11;
                dx(7) = w1*c32 - w3*c12;
                dx(8) = w1*c33 - w3*c13;
                dx(9) = w2*c11 - w1*c21;
                dx(10) = w2*c12 - w1*c22;
                dx(11) = w2*c13 - w1*c23;

                lielab::domain::halie _out{dx};
                return _out;
            };

            return out;
        }

        // vectorfield rigidbody(Eigen::MatrixXd inertia)
        // {
        //     // TODO: Rewrite this
        //     Eigen::MatrixXd invinertia = inertia.inverse();

        //     vectorfield out = [invinertia, inertia](const double _t, const lielab::domain::hmlie & _y)
        //     {
        //         lielab::domain::cola inp = std::get<lielab::domain::cola>(_y.space[0]);

        //         lielab::domain::so rbout(3);
        //         lielab::domain::cola force(3);
        //         lielab::domain::so dcmout(3);
        //         lielab::domain::su dq(2);

        //         double Ixx = inertia(0,0);

        //         rbout.set_vector(-invinertia*inp._data);

        //         double f = 0.0;
        //         if (_t >= 8.0 && _t <= 10.0)
        //         {
        //             f = 1.0;
        //         }
        //         force(0) = f;

        //         // Eigen::VectorXd temp(3);
        //         // temp(0) = inp._data(2)/2.0;
        //         // temp(1) = inp._data(1)/2.0;
        //         // temp(2) = inp._data(0)/2.0;

        //         dcmout.set_vector(inp._data);
        //         // dq.set_vector(temp);

        //         lielab::domain::halie out{rbout, dcmout};
        //         return out;
        //     };

        //     return out;
        // }

        Eigen::MatrixXd newton(Eigen::VectorXd x, Eigen::VectorXd F, double m, Eigen::VectorXd g)
        {
            /*!
            * Computes the rate-of-change using Newton's laws.
            */
            Eigen::VectorXd dx(3);

            dx(0) = F(0)/m - g(0);
            dx(1) = F(1)/m - g(1);
            dx(2) = F(2)/m - g(2);

            return dx;
        }
    }
}

#endif
