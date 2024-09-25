#ifndef _LIELAB_FUNCTIONS_COMMUTATOR_HPP
#define _LIELAB_FUNCTIONS_COMMUTATOR_HPP

#include "../abstract.hpp"
#include "../domain.hpp"

namespace Lielab
{
namespace functions
{
template <Lielab::abstract::LieAlgebra LA>
constexpr LA commutator(const LA & a, const LA & b)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    * This is the commutator function.
    */

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("commutator: Shapes of a and b must be equal.");
    }

    // Abelian speedhack. Return 0.
    if (a.abelian) return LA(a.shape);
    
    return a*b - b*a;
}
}
}

#endif
