#ifndef LIELAB_FUNCTIONS_LITTLECOAD_HPP_
#define LIELAB_FUNCTIONS_LITTLECOAD_HPP_

#include <vector>

#include "../abstract.hpp"
#include "../utils.hpp"
#include "commutator.hpp"
#include "littlead.hpp"

namespace Lielab
{
namespace functions
{

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::gl coad_numerical(const LA & a)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow \mathfrak{gl} \f}
    
    Numerical coadjoint function on Lie algebras.

    Arguments
    ---------
    @param[in] a First element.
    @param[out] out An instance of g.

    TODO
    ----
        - Compute wrt powers of coadjoint.
    
    */

    const Lielab::domain::gl ada = Lielab::functions::ad_numerical<LA>(a, 1);
    return Lielab::domain::gl(-ada.get_matrix().transpose());
}

template <Lielab::abstract::LieAlgebra LA>
LA coad_numerical(const LA & a, const LA & b)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    
    Numerical coadjoint function on Lie algebras.

    Arguments
    ---------
    @param[in] a First element.
    @param[out] out An instance of g.

    TODO
    ----
        - Compute wrt powers of coadjoint.
    
    */

    const Eigen::VectorXd bbar = b.get_vector();
    const Eigen::MatrixXd adsabar = Lielab::functions::coad_numerical<LA>(a).get_matrix();
    LA out = 0.0*b;
    out.set_vector(adsabar*bbar);
    return out;
}

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::gl coad(const LA & a)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow \mathfrak{gl} \f}
    
    Coadjoint function on Lie algebras.

    Arguments
    ---------
    @param[in] a First element.
    @param[out] out An instance of g.

    Notes
    -----

    Has some known shortcuts for performance.

    // 1. Power 0 is identity.

    // \f{equation*}{\text{ad}^{*0}_a = \mathbf{I}, \, \forall a \in \mathfrak{g} \f}

    2. Abelian Lie algebras return 0, except for power 0.

    \f{equation*}{\text{ad}^{*0}_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^{*j}_a = \mathbf{0}, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    TODO
    ----
        - This function generates the structure constants with each call, thus
          making it unsuitable for heavy numerical use.
        - There exist formula for adjoint power multiples of 2. This could
          accelerate this procedure by quite a bit at higher powers.
        - Compute wrt powers of coadjoint.
    
    */

    const ptrdiff_t dim = a.get_dimension();
    
    // Shortcut for power 0 adjoints.
    // if (p == 0)
    // {
    //     return Lielab::domain::gl(Eigen::MatrixXd::Identity(dim, dim));
    // }

    // Shortcut for Abelian Lie algebras.
    if (a.abelian)
    {
        return Lielab::domain::gl(dim);
    }

    const Lielab::domain::gl ada = Lielab::functions::ad<LA>(a, 1);
    return Lielab::domain::gl(-ada.get_matrix().transpose());
}

template <Lielab::abstract::LieAlgebra LA>
LA coad(const LA & a, const LA & b)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    
    Coadjoint function on Lie algebras.

    Arguments
    ---------
    @param[in] a First element.
    @param[out] out An instance of g.

    Notes
    -----

    Has some known shortcuts for performance.

    // 1. Power 0 is identity.

    // \f{equation*}{\text{ad}^{*0}_a = \mathbf{I}, \, \forall a \in \mathfrak{g} \f}

    2. Abelian Lie algebras return 0, except for power 0.

    \f{equation*}{\text{ad}^{*0}_a b = b \f}

    \f{equation*}{\text{ad}^{*j}_a b = b, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    TODO
    ----
        - This function generates the structure constants with each call, thus
          making it unsuitable for heavy numerical use.
        - There exist formula for adjoint power multiples of 2. This could
          accelerate this procedure by quite a bit at higher powers.
        - Compute wrt powers of coadjoint.
        - This uses g, not g* representation for argument 2. Whats the best way to handle this?
    
    */

    const ptrdiff_t dim = a.get_dimension();
    
    // Shortcut for power 0 adjoints.
    // if (p == 0)
    // {
    //     return Lielab::domain::gl(Eigen::MatrixXd::Identity(dim, dim));
    // }

    // Shortcut for Abelian Lie algebras.
    if (a.abelian)
    {
        return b;
    }

    const Eigen::VectorXd bbar = b.get_vector();
    const Eigen::MatrixXd adsabar = Lielab::functions::coad<LA>(a).get_matrix();
    LA out = 0.0*b;
    out.set_vector(adsabar*bbar);
    return out;
}

}
}

#endif
