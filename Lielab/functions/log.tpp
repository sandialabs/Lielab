#ifndef LIELAB_FUNCTIONS_LOG_TPP
#define LIELAB_FUNCTIONS_LOG_TPP

#include "log.hpp"

#include "Lielab/domain.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::functions
{

template<typename LG>
Lielab::domain::LieIII<LG> log_numerical(const LG& G)
{
    /*! \f{equation*}{ (G) \rightarrow \mathfrak{g} \f}
    
    This is the main logarithm function. Computes the logarithm in a
    computationally intensive numerical process.
    
    Arguments
    ---------
    @param[in] G An instance of G
    @param[out] out An instance of g

    */
    
    return (G.get_matrix()).log();
}

template<typename LG>
Lielab::domain::LieIII<LG> log(const LG& G)
{
    /*! \f{equation*}{ (G) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the exponential function. Will always use the
    numerical procedure.
    
    Arguments
    ---------
    @param[in] G An instance of G
    @param[out] out An instance of g

    */
    
    return log_numerical(G);
}

}

#endif
