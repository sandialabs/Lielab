#ifndef LIELAB_FUNCTIONS_EXPONENTIAL_TPP
#define LIELAB_FUNCTIONS_EXPONENTIAL_TPP

#include "exponential.hpp"

#include "adjoint.hpp"

#include "Lielab/domain.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::functions
{

template <typename g>
Lielab::domain::LieIII<g> exp_numerical(const g& x)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
    
    This is the main exponential function. Computes the exponential in a
    computationally intensive numerical process.
    
    Arguments
    ---------
    @param[in] x An instance of g
    @param[out] out An instance of G

    */

    using Lielab::domain::LieIII;
    return typename LieIII<g>::matrix_t((x.get_matrix()).exp());
}

template <typename g>
Lielab::domain::LieIII<g> exp(const g& x)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
    
    Catch-all function for the exponential function. Will always use the
    numerical procedure.
    
    Arguments
    ---------
    @param[in] x An instance of g
    @param[out] out An instance of G

    */
    
    return exp_numerical(x);
}

template <typename G>
Lielab::domain::LieIII<G> log_numerical(const G& g)
{
    /*! \f{equation*}{ (G) \rightarrow \mathfrak{g} \f}
    
    This is the main logarithm function. Computes the logarithm in a
    computationally intensive numerical process.
    
    Arguments
    ---------
    @param[in] g An instance of G
    @param[out] out An instance of g

    */
    
    using Lielab::domain::LieIII;
    return typename LieIII<G>::matrix_t((g.get_matrix()).log());
}

template <typename G>
Lielab::domain::LieIII<G> log(const G& g)
{
    /*! \f{equation*}{ (G) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the exponential function. Will always use the
    numerical procedure.
    
    Arguments
    ---------
    @param[in] g An instance of G
    @param[out] out An instance of g

    */
    
    return log_numerical(g);
}

template <typename g, typename out_t>
out_t dexp_numerical(const g& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    This is the main derivative of the exponential function. Computes it in the
    computationally intensive numerical process [1]:

    \f{equation*}{\text{dexp}_{x} = \sum_{j=0}^{\infty} \frac{1}{(j+1)!}\text{ad}_x^j \f}
    
    By default, we truncate at order 5 to best align with order 4 RK methods [2].

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    Notes
    -----

    Has a shortcut for performance.

    1. Abelian Lie algebras return identity for all orders.

    \f{equation*}{\text{dexp}_x = \sum_{j=0}^{\infty} \frac{1}{(j+1)!}\text{ad}^j_x = \text{ad}^0_x = \mathbf{I}, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    
    [2] Kenth Engø. On the construction of geometric integrators in the rkmk
        class. BIT Numerical Mathematics, 40:41–61, 2000.

    */

    out_t out = ad<g>(x, 0);

    // Special case where the domain is Abelian.
    if (x.is_abelian())
    {
        return out;
    }

    for (int ii = 1; ii <= order; ii++)
    {
        out += 1.0/std::tgamma(static_cast<double>(ii) + 2.0)*ad<g>(x, ii);
    }

    return out;
}

template <typename g, typename out_t>
out_t dexp(const g& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Catch-all function for the derivative of the exponential function.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a Instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexp_numerical<g>(x, order);
}

template <typename g>
g dexp_numerical(const g& x, const g& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    This is the main derivative of the exponential function. Computes it in the
    computationally intensive numerical process [1]:

    \f{equation*}{\text{dexp}_{x}(y) = \sum_{j=0}^{\infty} \frac{1}{(j+1)!}\text{ad}_x^j y \f}
    
    By default, we truncate at order 5 to best align with order 4 RK methods [2].

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] y Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    Notes
    -----

    Has a shortcut for performance.

    1. Abelian Lie algebras return b for all orders.

    \f{equation*}{\text{dexp}_x(y) = \sum_{j=0}^{\infty} \frac{1}{(j+1)!}\text{ad}^j_x(y) = \text{ad}^0_x(y) = y, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    
    [2] Kenth Engø. On the construction of geometric integrators in the rkmk
        class. BIT Numerical Mathematics, 40:41–61, 2000.

    */

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dexp_numerical: Shapes of a and b must be equal.");
    }

    g out(shape), adjc(shape);
    out = y;
    adjc = y;

    // Special case where the domain is abelian
    if (x.is_abelian())
    {
        return out;
    }

    for (int ii = 1; ii <= order; ii++)
    {
        adjc = commutator<g>(x, adjc);
        out += adjc*1.0/std::tgamma(static_cast<double>(ii) + 2.0);
    }

    return out;
}

template <typename g>
g dexp(const g& x, const g& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the derivative of the exponential function.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] y Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */
    
    return dexp_numerical<g>(x, y, order);
}

template <typename g, typename out_t>
out_t dlog_numerical(const g& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    This is the main derivative of the logarithm function. Computes
    it in a computationally intensive numerical process.

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance in the adjoint representation of g

    */

    return dexpinv_numerical<g>(x, order);
}

template <typename g, typename out_t>
out_t dlog(const g& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the derivative of the logarithm function.

    Arguments
    ---------
    @param[in] x Instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance in the adjoint representation of g

    */

    return dexpinv<g>(x, order);
}

template <typename g>
g dlog_numerical(const g& x, const g& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    This is the main derivative of the logarithm function. Computes
    it in a computationally intensive numerical process.

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] y Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexpinv_numerical<g>(x, y, order);
}

template <typename g>
g dlog(const g& x, const g& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the derivative of the logarithm function.

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] b Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexpinv<g>(x, y, order);
}

template <typename g, typename out_t>
out_t dloginv_numerical(const g& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    This is the main inverse derivative of the logarithm function. Computes
    it in a computationally intensive numerical process.

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance in the adjoint representation of g

    */

    return dexp_numerical<g>(x, order);
}

template <typename g, typename out_t>
out_t dloginv(const g& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the inverse derivative of the logarithm function.

    Arguments
    ---------
    @param[in] x Instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance in the adjoint representation of g

    */

    return dexp<g>(x, order);
}

template <typename g>
g dloginv_numerical(const g& x, const g& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    This is the main inverse derivative of the logarithm function. Computes
    it in a computationally intensive numerical process.

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] y Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexp_numerical<g>(x, y, order);
}

template <typename g>
g dloginv(const g& x, const g& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the inverse derivative of the logarithm function.

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] y Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexp<g>(x, y, order);
}

template <typename g, typename out_t>
out_t dexpinv_numerical(const g& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    This is the main inverse derivative of the exponential function. Computes
    it in the computationally intensive numerical process [1]:

    \f{equation*}{\text{dexp}_{a}^{-1} = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}^j_a \f}
    
    Where \f$B_j\f$ are Bernoulli numbers. By default, we truncate at order 5
    to best align with order 4 RK methods [2].

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance in the adjoint representation of g

    Notes
    -----

    Has a shortcut for performance.

    1. Abelian Lie algebras return identity for all orders.

    \f{equation*}{\text{dexpinv}_x = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}^j_x = \text{ad}^0_x = \mathbf{I}, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    
    [2] Kenth Engø. On the construction of geometric integrators in the rkmk
        class. BIT Numerical Mathematics, 40:41–61, 2000.

    */

    using Lielab::utils::bernoulli;

    out_t out = ad<g>(x, 0);

    // Special case where the domain is abelian
    if (x.is_abelian())
    {
        return out;
    }

    for (int ii = 1; ii <= order; ii++)
    {
        out += bernoulli(ii)/std::tgamma(static_cast<double>(ii) + 1.0)*ad<g>(x, ii);
    }

    return out;
}

template <typename g, typename out_t>
out_t dexpinv(const g& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the inverse derivative of the exponential function.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] x Instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance in the adjoint representation of g

    */

    return dexpinv_numerical<g>(x, order);
}

template <typename g>
g dexpinv_numerical(const g& x, const g& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    This is the main inverse derivative of the exponential function. Computes
    it in the computationally intensive numerical process [1]:

    \f{equation*}{\text{dexp}_{x}^{-1}(y) = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}_x^j y \f}
    
    Where \f$B_j\f$ are Bernoulli numbers. By default, we truncate at order 5
    to best align with order 4 RK methods [2].

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] y Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    Notes
    -----

    Has a shortcut for performance.

    1. Abelian Lie algebras return y for all orders.

    \f{equation*}{\text{dexpinv}_x(y) = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}^j_x(y) = \text{ad}^0_x(y) = y, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    
    [2] Kenth Engø. On the construction of geometric integrators in the rkmk
        class. BIT Numerical Mathematics, 40:41–61, 2000.

    */

    using Lielab::utils::bernoulli;

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dexpinv_numerical: Shapes of a and b must be equal.");
    }
    
    g out(shape), adjc(shape);
    adjc = y;
    out = adjc;

    // Special case where the domain is abelian
    if (x.is_abelian())
    {
        return out;
    }

    for (int ii = 1; ii <= order; ii++)
    {
        adjc = commutator<g>(x, adjc);
        out += adjc*bernoulli(ii)/std::tgamma(static_cast<double>(ii) + 1.0);
    }

    return out;
}

template <typename g>
g dexpinv(const g& x, const g& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the inverse derivative of the exponential function.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] x First instance of g
    @param[in] y Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexpinv_numerical<g>(x, y, order);
}

}

#endif
