#ifndef LIELAB_FUNCTIONS_LITTLECOAD_TPP
#define LIELAB_FUNCTIONS_LITTLECOAD_TPP

#include "commutator.hpp"
#include "littlead.hpp"

#include "Lielab/utils.hpp"

#include <vector>

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::glr coad_numerical(const LA & x, const int p)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow \mathfrak{glr} \f}
    
    Numerical coadjoint function on Lie algebras.

    Arguments
    ---------
    @param[in] x Lie algebra element.
    @param[in] p Power of the coadjoint.
    @param[out] out An instance of g.
    
    */

    const size_t dim = x.get_dimension();

    // Shortcut for power 0 coadjoints.
    if (p == 0)
    {
        return Lielab::domain::glr(Eigen::MatrixXd::Identity(dim, dim));
    }

    const Lielab::domain::glr ad1x = Lielab::functions::ad_numerical<LA>(x, 1);
    const Eigen::MatrixXd coad1xhat = -ad1x.get_matrix().transpose();

    // Raise coadjoint matrix to the specified power.
    Eigen::MatrixXd coadpxhat = Eigen::MatrixXd::Identity(dim, dim);
    for (size_t ii = 0; ii < p; ii++)
    {
        coadpxhat = coad1xhat*coadpxhat;
    }

    return Lielab::domain::glr(coadpxhat);
}

template <typename LA>
LA coad_numerical(const LA & x, const LA & y, const int p)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    
    Numerical coadjoint function on Lie algebras.

    Arguments
    ---------
    @param[in] x First element.
    @param[in] y Second element.
    @param[in] p Power of the coadjoint.
    @param[out] out An instance of g.
    
    */

    const Eigen::VectorXd ybar = y.get_vector();
    const Eigen::MatrixXd coadpxhat = Lielab::functions::coad_numerical<LA>(x, p).get_matrix();
    LA out = 0.0*y;
    out.set_vector(coadpxhat*ybar);
    return out;
}

template <typename LA>
Lielab::domain::glr coad(const LA & x, const int p)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow \mathfrak{glr} \f}
    
    Catch-all function for the coadjoint function on Lie algebras.
    Will sometimes use the numerical procedure.

    Arguments
    ---------
    @param[in] x First element.
    @param[in] p Power of the coadjoint.
    @param[out] out An instance of g.
    
    */

    const size_t dim = x.get_dimension();

    // Shortcut for power 0 coadjoints.
    if (p == 0)
    {
        return Lielab::domain::glr(Eigen::MatrixXd::Identity(dim, dim));
    }

    const Lielab::domain::glr ad1x = Lielab::functions::ad<LA>(x, 1);
    const Eigen::MatrixXd coad1xhat = -ad1x.get_matrix().transpose();

    // Raise coadjoint matrix to the specified power.
    Eigen::MatrixXd coadpxhat = Eigen::MatrixXd::Identity(dim, dim);
    for (size_t ii = 0; ii < static_cast<size_t>(p); ii++)
    {
        coadpxhat = coad1xhat*coadpxhat;
    }

    return Lielab::domain::glr(coadpxhat);
}

template <typename LA>
LA coad(const LA & x, const LA & y, const int p)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the coadjoint function on Lie algebras.
    Will sometimes use the numerical procedure.

    Arguments
    ---------
    @param[in] x First element.
    @param[in] y Second element.
    @param[out] out An instance of g.
    
    */

    const Eigen::VectorXd ybar = y.get_vector();
    const Eigen::MatrixXd coadpxhat = Lielab::functions::coad<LA>(x, p).get_matrix();
    LA out = 0.0*y;
    out.set_vector(coadpxhat*ybar);
    return out;
}

}

#endif
