#ifndef LIELAB_FUNCTIONS_AD_TPP
#define LIELAB_FUNCTIONS_AD_TPP

#include "littlead.hpp"

#include "exp.hpp"

#include "Lielab/domain.hpp"

#include <cassert>

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::GLR Ad_numerical(const LA & a)
{
    const Lielab::domain::glr ada = Lielab::functions::ad_numerical<LA>(a);
    return Lielab::functions::exp_numerical<Lielab::domain::glr>(ada);
}

template <typename LA>
Lielab::domain::GLR Ad(const LA & a)
{
    const Lielab::domain::glr ada = Lielab::functions::ad<LA>(a);
    return Lielab::functions::exp<Lielab::domain::glr>(ada);
}

// template <typename LG>
// Lielab::domain::GLR Ad_numerical(const LG & A)
// {
//     // TODO:
// }

// template <typename LG>
// Lielab::domain::GLR Ad(const LG & A)
// {
//     // TODO:
// }

template <typename LA>
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

    const size_t shape = a.get_shape();
    
    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("Ad: Shapes of a and b must be equal.");
    }

    const Lielab::domain::LieIII<LA> A = Lielab::functions::exp<LA>(a);
    return A.get_matrix()*b.get_matrix()*((A.inverse()).get_matrix());
}

template <typename LA>
LA Ad(const Lielab::domain::LieIII<LA> & A, const LA & b)
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
    
    const size_t shape = A.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("Ad: Shapes of A and b must be equal.");
    }
    
    return LA(A.get_matrix()*b.get_matrix()*((A.inverse()).get_matrix()));
}

}

#endif
