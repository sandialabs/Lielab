#ifndef _LIELAB_FUNCTIONS_DEXP_HPP
#define _LIELAB_FUNCTIONS_DEXP_HPP

namespace Lielab
{
    namespace functions
    {
        template <Lielab::abstract::LieAlgebra LA>
        LA dexp(const LA & a, const LA & b, const size_t order = 5)
        {
            /*!
            * Derivative of the exponential function.
            * 
            * Source: Iserles, A., Munthe-Kaas, H. Z., Nørsett, S. P., & Zanna, A. (2000).
            * Lie-group methods. Acta numerica, 9, 215-365.
            */

            LA out(b.shape), adjc(b.shape);
            out = b;
            adjc = b;

            // Special case where order = b returns b
            if (order == 0)
            {
                return out;
            }

            // Special case where the domain is abelian
            if (a.abelian)
            {
                return out;
            }

            for (int ii = 1; ii <= order; ii++)
            {
                adjc = commutator(a, adjc);
                out = out + adjc*1.0/factorial(ii+1);
                // TODO: The formula on wikipedia uses -1^ii but Munthe Kaas and Engo both don't have that term?
                // I don't know how to derive this myself and not sure of any identities to check.
                // out = out + adjc*std::pow(-1.0, static_cast<double>(ii))/factorial(ii+1);
            }

            // TODO: Munthe Kaas splits the G-action on g as d/dt exp(A) = dexp_(A)(d/dt A) exp(A)
            //       Do this operation inside or outside dexp? Need to solve the above TODO first though.
            // Lielab::domain::lieiii<LA> G = Lielab::functions::exp(a);
            // return LA(out.get_ados_representation()*G.get_ados_representation());

            return out;
        }
    }
}

#endif
