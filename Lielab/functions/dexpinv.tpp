#ifndef LIELAB_FUNCTIONS_DEXPINV_TPP
#define LIELAB_FUNCTIONS_DEXPINV_TPP

#include "dexpinv.hpp"

#include "bernoulli.hpp"
#include "dexp.hpp"
#include "littlead.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::glr dexpinv_numerical(const LA & a, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    This is the main inverse derivative of the exponential function. Computes
    it in the computationally intensive numerical process [1]:

    \f{equation*}{\text{dexp}_{a}^{-1} = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}^j_a \f}
    
    Where \f$B_j\f$ are Bernoulli numbers. By default, we truncate at order 5
    to best align with order 4 RK methods [2].

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    Notes
    -----

    Has a shortcut for performance.

    1. Abelian Lie algebras return identity for all orders.

    \f{equation*}{\text{dexpinv}_a = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}^j_a = \text{ad}^0_a = \mathbf{I}, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    
    [2] Kenth Engø. On the construction of geometric integrators in the rkmk
        class. BIT Numerical Mathematics, 40:41–61, 2000.

    */

    using Lielab::utils::factorial;

    Lielab::domain::glr out = ad<LA>(a, 0);

    // Special case where the domain is abelian
    if (a.abelian)
    {
        return out;
    }

    for (size_t ii = 1; ii <= order; ii++)
    {
        out = out + bernoulli(ii)/Lielab::utils::factorial(ii)*ad<LA>(a, ii);
    }

    return out;
}

template <typename LA>
Lielab::domain::glr dexpinv(const LA & a, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the inverse derivative of the exponential function.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a Instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexpinv_numerical<LA>(a, order);
}

template <typename LA>
LA dexpinv_numerical(const LA & a, const LA & b, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    This is the main inverse derivative of the exponential function. Computes
    it in the computationally intensive numerical process [1]:

    \f{equation*}{\text{dexp}_{a}^{-1}(b) = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}_a^j b \f}
    
    Where \f$B_j\f$ are Bernoulli numbers. By default, we truncate at order 5
    to best align with order 4 RK methods [2].

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] b Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    Notes
    -----

    Has a shortcut for performance.

    1. Abelian Lie algebras return b for all orders.

    \f{equation*}{\text{dexpinv}_a(b) = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}^j_a(b) = \text{ad}^0_a(b) = b, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    
    [2] Kenth Engø. On the construction of geometric integrators in the rkmk
        class. BIT Numerical Mathematics, 40:41–61, 2000.

    */

    using Lielab::utils::factorial;

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("dexpinv_numerical: Shapes of a and b must be equal.");
    }
    
    LA out(shape), adjc(shape);
    adjc = b;
    out = adjc;

    // Special case where the domain is abelian
    if (a.abelian)
    {
        return out;
    }

    for (size_t ii = 1; ii <= order; ii++)
    {
        adjc = commutator<LA>(a, adjc);
        out = out + adjc*bernoulli(ii)/factorial(ii);
    }

    return out;
}

template <typename LA>
LA dexpinv(const LA & a, const LA & b, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the inverse derivative of the exponential function.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] b Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexpinv_numerical<LA>(a, b, order);
}

}

#endif
