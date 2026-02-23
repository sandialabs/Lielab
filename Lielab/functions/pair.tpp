#ifndef LIELAB_FUNCTIONS_PAIR_TPP
#define LIELAB_FUNCTIONS_PAIR_TPP

#include "pair.hpp"

//Z#include "Lielab/testing.hpp"

namespace Lielab::functions
{

template <typename LA>
double pair(const LA& a, const LA& b)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathbb{R} \f}
    * 
    * Pairing of vectors on themself.
    *
    * TODO: Check shapes or dimension?
    */


    //Zlielab_assert(a.get_shape() == b.get_shape(), "Shapes must be equal.");
    return a.get_vector().dot(b.get_vector());
}

}

#endif
