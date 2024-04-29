#ifndef _LIELAB_FUNCTIONS_DEXPINV_HPP
#define _LIELAB_FUNCTIONS_DEXPINV_HPP

namespace Lielab
{
    namespace functions
    {
        template <Lielab::abstract::LieAlgebra LA>
        LA dexpinv(const LA & a, const LA & b, const size_t order = 5)
        {
            /*!
            * Derivative of the inverse exponential function.
            * 
            * Source: Iserles, A., Munthe-Kaas, H. Z., Nørsett, S. P., & Zanna, A. (2000).
            * Lie-group methods. Acta numerica, 9, 215-365.
            */
            
            LA out(b.shape), adjc(b.shape);

            // Special case where order = 0 returns 0 vector.
            if (order == 0)
            {
                return out;
            }

            adjc = b;
            out = adjc;

            // Special case where the domain is abelian
            if (a.abelian)
            {
                return out;
            }

            for (int ii = 1; ii < order; ii++)
            {
                adjc = commutator(a, adjc);
                out = out + adjc*bernoulli(ii)/factorial(ii);
            }
            return out;
        }
    }
}

#endif
