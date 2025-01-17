#ifndef _LIELAB_FUNCTIONS_DLOG_HPP
#define _LIELAB_FUNCTIONS_DLOG_HPP

#include "dexpinv.hpp"

namespace Lielab
{
namespace functions
{

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::gl dlog_numerical(const LA & a, const size_t order = 5)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
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

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::gl dlog(const LA & a, const size_t order = 5)
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

template <Lielab::abstract::LieAlgebra LA>
LA dlog_numerical(const LA & a, const LA & b, const size_t order = 5)
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

template <Lielab::abstract::LieAlgebra LA>
LA dlog(const LA & a, const LA & b, const size_t order = 5)
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
}

#endif
