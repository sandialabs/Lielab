#ifndef LIELAB_FUNCTIONS_ADJOINT_TPP
#define LIELAB_FUNCTIONS_ADJOINT_TPP

#include "adjoint.hpp"

#include "Lielab/utils.hpp"

#include <Eigen/Core>

#include <vector>

namespace Lielab::functions
{

template <typename g>
g commutator(const g& a, const g& b)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    * This is the commutator function.
    */

    const int shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("commutator: Shapes of a and b must be equal.");
    }

    // Abelian speedhack. Return 0.
    if (a.is_abelian())
    {
        return g::zero(shape);
    }

    const typename g::matrix_t ahat = a.get_matrix();
    const typename g::matrix_t bhat = b.get_matrix();
    
    return g(ahat*bhat - bhat*ahat);
}

template <typename g, typename out_t>
out_t ad_numerical(const g& a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Adjoint function on Lie algebras.
    
    This computes \f$\text{ad}_a^p\f$ in an entirely numerical manner. For basis
    \f$e_i \in \mathfrak{g}\f$, the Lie algebra structure constants are generated with:

    \f{equation*}{[e_i, e_j] = C_{ijk} e_k\f}

    Then the first adjoint of \f$a\f$ is found with:

    \f{equation*}{[\text{ad}_a]_{kj} = C_{ajk} \f}

    Powers of the adjoint are then simply

    \f{equation*}{\text{ad}^p_a = \Pi_p \text{ad}_a \f}

    Arguments
    ---------
    @param[in] a First element.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of g.

    Notes
    -----

    Has some known shortcuts for performance.

    1. Power 0 is identity.

    \f{equation*}{\text{ad}^0_a = \mathbf{I}, \, \forall a \in \mathfrak{g} \f}

    2. Abelian Lie algebras return 0, except for power 0.

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^j_a = \mathbf{0}, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    TODO
    ----
        - This function generates the structure constants with each call, thus
          making it unsuitable for heavy numerical use.
        - There exist formula for adjoint power multiples of 2. This could
          accelerate this procedure by quite a bit at higher powers.
    
    */

    const int dim = a.get_dimension();
    
    // Shortcut for power 0 adjoints.
    if (p == 0)
    {
        return out_t(Eigen::MatrixXd::Identity(dim, dim));
    }

    // Shortcut for Abelian Lie algebras.
    if (a.is_abelian())
    {
        return out_t(Eigen::MatrixXd::Zero(dim, dim));
    }

    const Eigen::VectorXd abar = a.get_vector();
    const int shape = a.get_shape();
    
    // Generate basis.
    std::vector<g> basis;
    for (int ii = 0; ii < dim; ii++)
    {
        basis.push_back(g::basis(ii, shape));
    }

    // Generate structure constants.
    std::vector<std::vector<Eigen::VectorXd>> C;
    for (int ii = 0; ii < dim; ii++)
    {
        std::vector<Eigen::VectorXd> vecj;
        for (int jj = 0; jj < dim; jj++)
        {
            const Eigen::VectorXd veck = commutator(basis[ii], basis[jj]).get_vector();
            vecj.push_back(veck);
        }
        C.push_back(vecj);
    }

    // Generate first adjoint matrix.
    Eigen::MatrixXd adja = Eigen::MatrixXd::Zero(dim, dim);
    for (int ii = 0; ii < dim; ii++)
    {
        for (int jj = 0; jj < dim; jj++)
        {
            for (int kk = 0; kk < dim; kk++)
            {
                adja(jj, ii) += abar(kk)*C[kk][ii](jj);
            }
        }
    }

    // Raise adjoint matrix to the specified power.
    const Eigen::MatrixXd adjap = adja.pow(p);

    return out_t(adjap);
}

template <typename g, typename out_t>
out_t ad(const g& a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Catch-all function for the adjoint function on Lie algebras.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a First element.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of g.

    */
    
    return ad_numerical<g>(a, p);
}

template <typename g>
g ad_numerical(const g& a, const g& b, const int p)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Adjoint function on Lie algebras.
    
    This computes \f$\text{ad}_a^p b\f$ in an entirely numerical manner:

    \f{equation*}{\text{ad}^0_a b = b \f}

    \f{equation*}{\text{ad}^1_a b = [a, b] \f}

    \f{equation*}{\text{ad}^2_a b = [a, [a, b]] \f}

    ... and so on.

    Arguments
    ---------
    @param[in] a First element.
    @param[in] b Second element.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of g.
    
    */

    const int shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("ad_numerical: Shapes of a and b must be equal.");
    }

    g adjb = b;

    for (ptrdiff_t ii = 0; ii < p; ii++)
    {
        adjb = commutator(a, adjb);
    }

    return adjb;
}

template <typename g>
g ad(const g& a, const g& b, const int p)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the adjoint function on Lie algebras.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a First element.
    @param[in] b Second element.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of g.

    */
    
    return ad_numerical<g>(a, b, p);
}

// template <typename LA>
// Lielab::domain::GLR Ad_numerical(const LA& a)
// {
//     const Lielab::domain::glr ada = Lielab::functions::ad_numerical<LA>(a);
//     return Lielab::functions::exp_numerical<Lielab::domain::glr>(ada);
// }

// template <typename LA>
// Lielab::domain::GLR Ad(const LA& a)
// {
//     const Lielab::domain::glr ada = Lielab::functions::ad<LA>(a);
//     return Lielab::functions::exp<Lielab::domain::glr>(ada);
// }

// template <typename LG>
// Lielab::domain::GLR Ad_numerical(const LG& A)
// {
//     // TODO:
// }

// template <typename LG>
// Lielab::domain::GLR Ad(const LG& A)
// {
//     // TODO:
// }

// template <typename LA>
// LA Ad(const LA& a, const LA& b)
// {
//     /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
//     * 
//     * Group Adjoint action.
//     * 
//     * \f{equation*}{ \text{Ad}(a, b) = exp(a) b exp(a)^{-1} \f}
//     * 
//     * @param[in] a A Lie algebra.
//     * @param[in] b A Lie algebra.
//     * @param[out] out A Lie algebra.
//     */

//     const size_t shape = a.get_shape();
    
//     if (shape != b.get_shape())
//     {
//         throw Lielab::utils::InputError("Ad: Shapes of a and b must be equal.");
//     }

//     const Lielab::domain::LieIII<LA> A = Lielab::functions::exp<LA>(a);
//     return A.get_matrix()*b.get_matrix()*((A.inverse()).get_matrix());
// }

template <typename G>
Lielab::domain::LieIII<G> Ad(const G& g, const Lielab::domain::LieIII<G>& x)
{
    /*! \f{equation*}{ (G, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    * 
    * Group Adjoint action.
    * 
    * \f{equation*}{ \text{Ad}(g, x) = gxg^{-1} \f}
    * 
    * @param[in] g A Lie group.
    * @param[in] x A Lie algebra.
    * @param[out] out A Lie algebra.
    */
    
    const int shape = g.get_shape();

    if (shape != x.get_shape())
    {
        throw Lielab::utils::InputError("Ad: Shapes of g and x must be equal.");
    }
    
    return Lielab::domain::LieIII<G>(g.get_matrix()*x.get_matrix()*((g.inverse()).get_matrix()));
}

template <typename g, typename out_t>
out_t coad_numerical(const g& x, const int power)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow \mathfrak{glr} \f}
    
    Numerical coadjoint function on Lie algebras.

    Arguments
    ---------
    @param[in] x Lie algebra element.
    @param[in] power Power of the coadjoint.
    @param[out] out An instance in the coadjoint representation of g.
    
    */

    const size_t dim = x.get_dimension();

    // Shortcut for power 0 coadjoints.
    if (power == 0)
    {
        return out_t(Eigen::MatrixXd::Identity(dim, dim));
    }

    const out_t ad1x = Lielab::functions::ad_numerical<g>(x, 1);
    const Eigen::MatrixXd coad1xhat = -ad1x.get_matrix().transpose();

    // Raise coadjoint matrix to the specified power.
    Eigen::MatrixXd coadpxhat = Eigen::MatrixXd::Identity(dim, dim);
    for (int ii = 0; ii < power; ii++)
    {
        coadpxhat = coad1xhat*coadpxhat;
    }

    return out_t(coadpxhat);
}

template <typename g>
g coad_numerical(const g& x, const g& y, const int power)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    
    Numerical coadjoint function on Lie algebras.

    Arguments
    ---------
    @param[in] x First element.
    @param[in] y Second element.
    @param[in] power Power of the coadjoint.
    @param[out] out An instance of g.
    
    */

    const Eigen::VectorXd ybar = y.get_vector();
    const Eigen::MatrixXd coadpxhat = Lielab::functions::coad_numerical<g>(x, power).get_matrix();
    g out = 0.0*y;
    out.set_vector(coadpxhat*ybar);
    return out;
}

template <typename g, typename out_t>
out_t coad(const g& x, const int power)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow \mathfrak{glr} \f}
    
    Catch-all function for the coadjoint function on Lie algebras.
    Will sometimes use the numerical procedure.

    Arguments
    ---------
    @param[in] x First element.
    @param[in] p Power of the coadjoint.
    @param[out] out An instance in the coadjoint representation of g.
    
    */

    const size_t dim = x.get_dimension();

    // Shortcut for power 0 coadjoints.
    if (power == 0)
    {
        return out_t(Eigen::MatrixXd::Identity(dim, dim));
    }

    const out_t ad1x = Lielab::functions::ad<g>(x, 1);
    const Eigen::MatrixXd coad1xhat = -ad1x.get_matrix().transpose();

    // Raise coadjoint matrix to the specified power.
    Eigen::MatrixXd coadpxhat = Eigen::MatrixXd::Identity(dim, dim);
    for (int ii = 0; ii < power; ii++)
    {
        coadpxhat = coad1xhat*coadpxhat;
    }

    return out_t(coadpxhat);
}

template <typename g>
g coad(const g& x, const g& y, const int power)
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
    const Eigen::MatrixXd coadpxhat = Lielab::functions::coad<g>(x, power).get_matrix();
    g out = 0.0*y;
    out.set_vector(coadpxhat*ybar);
    return out;
}

}

#endif
