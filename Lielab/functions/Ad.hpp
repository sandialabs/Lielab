#ifndef _LIELAB_FUNCTIONS_AD_HPP
#define _LIELAB_FUNCTIONS_AD_HPP

#include "../abstract.hpp"
#include "../domain.hpp"

#include <cassert>

namespace Lielab
{
    namespace functions
    {
        template <Lielab::abstract::LieAlgebra LA>
        LA Ad(const Lielab::domain::lieiii<LA> & G, const LA & g)
        {
            /*! \f{equation*}{ (G, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
            * 
            * Group Adjoint action.
            * 
            * \f{equation*}{ Ad(G, X) = GXG^{-1} \f}
            * 
            * @param[in] G A Lie group.
            * @param[in] g A lie algebra.
            * @param[out] out A lie algebra.
            */
            
            assert(G.shape == g.shape);
            return G.get_matrix()*g.get_matrix()*(G.get_matrix().inverse());
        }
    }
}

#endif
