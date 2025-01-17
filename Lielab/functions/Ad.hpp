#ifndef _LIELAB_FUNCTIONS_AD_HPP
#define _LIELAB_FUNCTIONS_AD_HPP

#include "../abstract.hpp"
#include "../domain.hpp"

#include "littlead.hpp"
#include "exp.hpp"

#include <cassert>

namespace Lielab
{
namespace functions
{

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::GL Ad_numerical(const LA & a)
{
    const Lielab::domain::gl ada = Lielab::functions::ad_numerical<LA>(a);
    return Lielab::functions::exp_numerical<Lielab::domain::gl>(ada);
}

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::GL Ad(const LA & a)
{
    const Lielab::domain::gl ada = Lielab::functions::ad<LA>(a);
    return Lielab::functions::exp<Lielab::domain::gl>(ada);
}

// template <Lielab::abstract::LieGroup LG>
// Lielab::domain::GL Ad_numerical(const LG & A)
// {
//     // TODO:
// }

// template <Lielab::abstract::LieGroup LG>
// Lielab::domain::GL Ad(const LG & A)
// {
//     // TODO:
// }

template <Lielab::abstract::LieAlgebra LA>
LA Ad(const LA & a, const LA & b)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    * 
    * Group Adjoint action.
    * 
    * \f{equation*}{ \text{Ad}(a, b) = exp(a) b exp(a)^{-1} \f}
    * 
    * @param[in] a A Lie algebra.
    * @param[in] b A Lie algebra.
    * @param[out] out A Lie algebra.
    */
    
    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("Ad: Shapes of a and b must be equal.");
    }

    const Lielab::domain::lieiii<LA> A = Lielab::functions::exp<LA>(a);
    return A.get_matrix()*b.get_matrix()*((A.inverse()).get_matrix());
}

template <Lielab::abstract::LieAlgebra LA>
LA Ad(const Lielab::domain::lieiii<LA> & A, const LA & b)
{
    /*! \f{equation*}{ (G, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    * 
    * Group Adjoint action.
    * 
    * \f{equation*}{ \text{Ad}(A, b) = AbA^{-1} \f}
    * 
    * @param[in] A A Lie group.
    * @param[in] b A Lie algebra.
    * @param[out] out A Lie algebra.
    */
    
    if (A.shape != b.shape)
    {
        throw Lielab::utils::InputError("Ad: Shapes of A and b must be equal.");
    }
    
    return A.get_matrix()*b.get_matrix()*((A.inverse()).get_matrix());
}
}
}

#endif
