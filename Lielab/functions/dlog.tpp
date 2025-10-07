#ifndef LIELAB_FUNCTIONS_DLOG_TPP
#define LIELAB_FUNCTIONS_DLOG_TPP

#include "dlog.hpp"

#include "dexpinv.hpp"

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::glr dlog_numerical(const LA & a, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    This is the main derivative of the logarithm function. Computes
    it in a computationally intensive numerical process.

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexpinv_numerical<LA>(a, order);
}

template <typename LA>
Lielab::domain::glr dlog(const LA & a, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the derivative of the logarithm function.

    Arguments
    ---------
    @param[in] a Instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexpinv<LA>(a, order);
}

template <typename LA>
LA dlog_numerical(const LA & a, const LA & b, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    This is the main derivative of the logarithm function. Computes
    it in a computationally intensive numerical process.

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] b Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexpinv_numerical<LA>(a, b, order);
}

template <typename LA>
LA dlog(const LA & a, const LA & b, const size_t order)
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

    return dexpinv<LA>(a, b, order);
}

}

#endif
