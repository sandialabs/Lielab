#ifndef _LIELAB_FUNCTIONS_H
#define _LIELAB_FUNCTIONS_H

#include <stdexcept>
#include <iostream>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <complex>

#include "abstract.hpp"
#include "domain"

using std::sin;
using std::cos;
using std::pow;

using lielab::domain::lieiii;

namespace lielab
{
    /*!
    * 
    */
    namespace functions
    {
        template<lielab::abstract::LieAlgebra LA>
        LA ad(const LA & a, const LA & b)
        {
            /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
            * 
            * adjoint function.
            */

            return commutator(a, b);
        }

        // template<lielab::abstract::LieAlgebra LA>
        // lielab::domain::cola coad(const LA & g, const lielab::domain::cola & mu)
        // {
        //     /*! \f{equation}{ (\mathfrak{g}, \mathfrak{g}^*) \rightarrow \mathfrak{g}^* }
        //     *
        //     * Coadjoint function.
        //     */

        //     std::vector<LA> basis;
        //     const size_t d1 = g.get_dimension();
        //     const size_t d2 = mu.get_dimension();
        //     // double norm = 0.0;

        //     if (d1 != d2)
        //     {
        //         throw SizeError("Input data malformed: dimension " + std::to_string(d1) + " != " + std::to_string(d2) + ".");
        //     }
            
        //     Eigen::VectorXd out(d1);

        //     const size_t s = g.shape;

        //     for (int ii = 0; ii < d1; ii++)
        //     {
        //         basis.push_back(LA::basis(ii, s));
        //     }

        //     for (int ii = 0; ii < d1; ii++)
        //     {
        //         out(ii) = copair(mu, -ad(g, basis[ii]));
        //     }

        //     return out;
        // }
    }
}

#endif
