#ifndef LIELAB_FUNCTIONS_PAIR_TPP
#define LIELAB_FUNCTIONS_PAIR_TPP

#include "pair.hpp"

#include <cassert>

namespace Lielab::functions
{

template <typename LA>
double pair(const LA & a, const LA & b)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathbb{R} \f}
    * 
    * Pairing of vectors on themself.
    *
    * TODO: Check shapes or dimension?
    */

    if (a.get_dimension() != b.get_dimension())
    {
        throw Lielab::utils::InputError("pair: Dimensions of a and b must be equal.");
    }

    return a.get_vector().dot(b.get_vector());
}

}

#endif
