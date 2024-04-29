#ifndef _LIELAB_FUNCTIONS_EXP_HPP
#define _LIELAB_FUNCTIONS_EXP_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace functions
    {
        template<typename LA>
        Lielab::domain::lieiii<LA> _exp(const LA & la)
        {
            /*! \f{equation*}{ (\mathfrak{g} \rightarrow G \f}
            *
            * This is the main exponential function. Computes
            * the exponential in a computationally intensive
            * numerical process.
            */

            return (la.get_matrix()).exp();
        }

        template<typename LA>
        Lielab::domain::lieiii<LA> exp(const LA & la)
        {
            /*! \f{equation*}{ (\mathfrak{g} \rightarrow G \f}
            *
            * This is the exponential function.
            */
            
            return _exp(la);
        }

        template<>
        Lielab::domain::CN exp(const Lielab::domain::cn & la)
        {
            /*! \f{equation*}{ \mathfrak{cn} \rightarrow CN \f}
            *
            * cn exponential overload
            */

            return la._data;
        }

        template<>
        Lielab::domain::RN exp(const Lielab::domain::rn & la)
        {
            /*! \f{equation*}{ \mathfrak{rn} \rightarrow RN \f}
            *
            * rn exponential overload
            */

            return la._data;
        }

        template<>
        Lielab::domain::SO exp(const Lielab::domain::so & la)
        {
            /*! \f{equation*}{ \mathfrak{so} \rightarrow SO \f}
            *
            * so exponential overload
            *
            * Uses the Euler-Rodrigues formula for so(3), otherwise
            * calls the numerical method for so(n), n != 3.
            */

            static const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(3,3);
            static const Eigen::MatrixXd ZZ = Eigen::MatrixXd::Zero(3,3);

            if (la.shape == 3)
            {
                Lielab::domain::SO out(3);

                const double beta2 = la(0,1)*la(0,1) + la(0,2)*la(0,2) + la(1,2)*la(1,2);
                const double beta = std::sqrt(beta2);

                if (beta < 1e-14)
                {
                    return out;
                }

                Eigen::MatrixXd S1 = Id;
                Eigen::MatrixXd S2 = ZZ;

                const double a = std::sin(beta)/beta;
                const double b = std::sin(beta/2);
                const double c = 2*b*b/beta2;

                S1(0,1) = a*la(0,1);
                S1(0,2) = a*la(0,2);
                S1(1,2) = a*la(1,2);
                S1(1,0) = -S1(0,1);
                S1(2,0) = -S1(0,2);
                S1(2,1) = -S1(1,2);
                
                const double cS12 = c*la(0,1);
                const double cS13 = c*la(0,2);
                const double v2 = cS12*la(0,1);
                const double u2 = cS13*la(0,2);
                const double t2 = c*la(1,2)*la(1,2);

                S2(0,0) = -(v2+u2);
                S2(0,1) = -cS13*la(1,2);
                S2(0,2) = cS12*la(1,2);
                S2(1,0) = S2(0,1);
                S2(1,1) = -(v2+t2);
                S2(1,2) = -cS12*la(0,2);
                S2(2,0) = S2(0,2);
                S2(2,1) = S2(1,2);
                S2(2,2) = -(u2+t2);

                out._data = S1+S2;
                return out;
            }

            return _exp(la);
        }
    }
}

#endif
