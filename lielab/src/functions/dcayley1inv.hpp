#ifndef _LIELAB_FUNCTIONS_DCAYLEY1INV_HPP
#define _LIELAB_FUNCTIONS_DCAYLEY1INV_HPP

namespace lielab
{
    namespace functions
    {
        template <lielab::abstract::LieAlgebra LA>
        LA dcayley1inv(const LA & u, const LA & v)
        {
            /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
            * 
            * Derivative of the Cayley function.
            *
            * \f{equation*}{ \text{dcay}_u^{-1}(v) = v - \frac{1}{2} [u,v] - \frac{1}{4} u \cdot v \cdot u \f}
            * 
            * Source: Engø, Kenth. "On the construction of geometric integrators in the RKMK class."
            * BIT Numerical Mathematics 40.1 (2000): 41-61.
            */

            LA temp(v.shape);
            temp = v;

            if (!v.abelian)
            {
                temp = temp - 1.0/2.0*commutator(u, v);
            }

            return temp - 1.0/4.0*u*v*u;
        }
    }
}

#endif
