#ifndef _LIELAB_RK_METHODS_HPP
#define _LIELAB_RK_METHODS_HPP


#include <array>

#include <Eigen/Core>

namespace Lielab
{
    namespace topos
    {
        enum class RKMETHOD {E1, RK45};

        Eigen::MatrixXd RKMETHOD_to_A(const RKMETHOD & rkmethod)
        {
            if (rkmethod == RKMETHOD::E1)
            {
                Eigen::MatrixXd A(1,1);
                A(0,0) = 0.0;
                return A;
            }

            if (rkmethod == RKMETHOD::RK45)
            {
                Eigen::MatrixXd A(6,6);
                A <<    0.0,            0.0,           0.0,           0.0,         0.0,       0.0,
                        1.0/4.0,        0.0,           0.0,           0.0,         0.0,       0.0,
                        3.0/32.0,       9.0/32.0,      0.0,           0.0,         0.0,       0.0,
                     1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0,    0.0,         0.0,       0.0,
                      439.0/216.0,     -8.0,        3680.0/513.0,  -845.0/4104.0,  0.0,       0.0,
                       -8.0/27.0,       2.0,       -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0, 0.0;
                return A;
            }

            return Eigen::MatrixXd::Zero(0,0);
        }

        Eigen::VectorXd RKMETHOD_to_B(const RKMETHOD & rkmethod)
        {
            if (rkmethod == RKMETHOD::E1)
            {
                Eigen::VectorXd B(1);
                B(0) = 1.0;
                return B;
            }

            if (rkmethod == RKMETHOD::RK45)
            {
                Eigen::VectorXd B(6);
                B << 25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0;
                return B;
            }

            return Eigen::VectorXd::Zero(0);
        }

        Eigen::VectorXd RKMETHOD_to_Bhat(const RKMETHOD & rkmethod)
        {
            if (rkmethod == RKMETHOD::E1)
            {
                Eigen::VectorXd Bhat(1);
                Bhat(0) = 0.0;
                return Bhat;
            }

            if (rkmethod == RKMETHOD::RK45)
            {
                Eigen::VectorXd Bhat(6);
                Bhat << 16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0;
                return Bhat;
            }

            return Eigen::VectorXd::Zero(0);
        }

        Eigen::VectorXd RKMETHOD_to_C(const RKMETHOD & rkmethod)
        {
            if (rkmethod == RKMETHOD::E1)
            {
                Eigen::VectorXd C(1);
                C(0) = 0.0;
                return C;
            }

            if (rkmethod == RKMETHOD::RK45)
            {
                Eigen::VectorXd C(6);
                C << 0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0;
                return C;
            }

            return Eigen::VectorXd::Zero(0);
        }

        constexpr size_t RKMETHOD_to_order(const RKMETHOD & rkmethod)
        {
            if (rkmethod == RKMETHOD::E1) return 1;
            if (rkmethod == RKMETHOD::RK45) return 4;

            return 0;
        }

        constexpr size_t RKMETHOD_to_n(const RKMETHOD & rkmethod)
        {
            if (rkmethod == RKMETHOD::E1) return 1;
            if (rkmethod == RKMETHOD::RK45) return 6;

            return 0;
        }

        constexpr bool RKMETHOD_to_variable(const RKMETHOD & rkmethod)
        {
            if (rkmethod == RKMETHOD::E1) return false;
            if (rkmethod == RKMETHOD::RK45) return true;

            return false;
        }
    }
}

#endif
