#ifndef LIELAB_FUNCTIONS_COMMUTATOR_TPP
#define LIELAB_FUNCTIONS_COMMUTATOR_TPP

#include "commutator.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <typename LA>
LA commutator(const LA& a, const LA& b)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    * This is the commutator function.
    */

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("commutator: Shapes of a and b must be equal.");
    }

    // Abelian speedhack. Return 0.
    if (a.abelian) return LA(shape); // TODO: Make this ::from_shape

    const typename LA::matrix_t ahat = a.get_matrix();
    const typename LA::matrix_t bhat = b.get_matrix();
    
    return LA(ahat*bhat - bhat*ahat);
}

}

#endif
