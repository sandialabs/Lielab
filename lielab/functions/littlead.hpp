#ifndef LIELAB_FUNCTIONS_LITTLEAD_HPP_
#define LIELAB_FUNCTIONS_LITTLEAD_HPP_

#include "../abstract/abstract_all.hpp"
#include "commutator.hpp"

namespace Lielab
{
/*!
* 
*/
namespace functions
{

template<Lielab::abstract::LieAlgebra LA>
LA ad(const LA & a, const LA & b)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    * 
    * (little) adjoint function.
    */

    return commutator(a, b);
}

}
}

#endif
