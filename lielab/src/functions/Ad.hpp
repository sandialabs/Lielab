#ifndef _LIELAB_FUNCTIONS_AD_HPP
#define _LIELAB_FUNCTIONS_AD_HPP

#include <cassert>

namespace lielab
{
    namespace functions
    {
        template <lielab::abstract::LieAlgebra LA>
        LA Ad(const lielab::domain::lieiii<LA> & G, const LA & g)
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
            return G.get_ados_representation()*g.get_ados_representation()*G.get_ados_representation().inverse();
        }
    }
}

#endif
