#ifndef _LIELAB_FUNCTIONS_PAIR_HPP
#define _LIELAB_FUNCTIONS_PAIR_HPP

#include <cassert>

namespace lielab
{
    namespace functions
    {
        template <lielab::abstract::LieAlgebra LA>
        double pair(const LA & a, const LA & b)
        {
            /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathbb{R} \f}
            * 
            * Pairing of vectors on themself.
            */

            assert(a.get_dimension() == b.get_dimension());
            return a.get_vector().dot(b.get_vector());
        }
    }
}

#endif
