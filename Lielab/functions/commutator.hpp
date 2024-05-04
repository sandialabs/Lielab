#ifndef _LIELAB_FUNCTIONS_COMMUTATOR_HPP
#define _LIELAB_FUNCTIONS_COMMUTATOR_HPP

#include "../abstract.hpp"
#include "../domain.hpp"

namespace Lielab
{
    namespace functions
    {
        template <Lielab::abstract::LieAlgebra LA>
        constexpr LA commutator(const LA & g1, const LA & g2)
        {
            /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
            * This is the commutator function.
            */

            // Abelian speedhack. Return 0.
            if (g1.abelian) return LA(g1.shape);
            
            return g1*g2 - g2*g1;
        }
    }
}

#endif
