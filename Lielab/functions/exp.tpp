#ifndef LIELAB_FUNCTIONS_EXP_TPP
#define LIELAB_FUNCTIONS_EXP_TPP

#include "exp.hpp"

#include "Lielab/domain.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::functions
{

template<typename LA>
Lielab::domain::LieIII<LA> exp_numerical(const LA & la)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
    
    This is the main exponential function. Computes the exponential in a
    computationally intensive numerical process.
    
    Arguments
    ---------
    @param[in] la An instance of g
    @param[out] out An instance of G

    */

    return (la.get_matrix()).exp();
}

template<typename LA>
Lielab::domain::LieIII<LA> exp(const LA & la)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
    
    Catch-all function for the exponential function. Will always use the
    numerical procedure.
    
    Arguments
    ---------
    @param[in] la An instance of g
    @param[out] out An instance of G

    */
    
    return exp_numerical(la);
}

}

#endif
