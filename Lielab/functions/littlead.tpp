#ifndef LIELAB_FUNCTIONS_LITTLEAD_TPP
#define LIELAB_FUNCTIONS_LITTLEAD_TPP

#include "littlead.hpp"

#include "commutator.hpp"

#include "Lielab/utils.hpp"

#include <Eigen/Core>

#include <vector>

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::glr ad_numerical(const LA & a, const int p)
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

    const size_t dim = a.get_dimension();
    
    // Shortcut for power 0 adjoints.
    if (p == 0)
    {
        return Lielab::domain::glr(Eigen::MatrixXd::Identity(dim, dim));
    }

    // Shortcut for Abelian Lie algebras.
    if (a.abelian)
    {
        return Lielab::domain::glr(dim);
    }

    const Eigen::VectorXd ahat = a.get_vector();
    const size_t shape = a.get_shape();
    
    // Generate basis.
    std::vector<LA> basis;
    for (size_t ii = 0; ii < dim; ii++)
    {
        basis.push_back(LA::basis(ii, shape));
    }

    // Generate structure constants.
    std::vector<std::vector<Eigen::VectorXd>> C;
    for (size_t ii = 0; ii < dim; ii++)
    {
        std::vector<Eigen::VectorXd> vecj;
        for (size_t jj = 0; jj < dim; jj++)
        {
            const Eigen::VectorXd veck = commutator(basis[ii], basis[jj]).get_vector();
            vecj.push_back(veck);
        }
        C.push_back(vecj);
    }

    // Generate first adjoint matrix.
    Eigen::MatrixXd adja = Eigen::MatrixXd::Zero(dim, dim);
    for (size_t ii = 0; ii < dim; ii++)
    {
        for (size_t jj = 0; jj < dim; jj++)
        {
            for (size_t kk = 0; kk < dim; kk++)
            {
                adja(jj, ii) += ahat(kk)*C[kk][ii](jj);
            }
        }
    }

    // Raise adjoint matrix to the specified power.
    // Eigen::MatrixXd adjap = Eigen::MatrixXd::Identity(dim, dim);
    // for (size_t ii = 0; ii < static_cast<size_t>(p); ii++)
    // {
    //     adjap = adja*adjap;
    // }
    const Eigen::MatrixXd adjap = adja.pow(p);

    return Lielab::domain::glr(adjap);
}

template <typename LA>
Lielab::domain::glr ad(const LA & a, const int p)
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
    
    return ad_numerical<LA>(a, p);
}

template <typename LA>
LA ad_numerical(const LA & a, const LA & b, const int p)
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

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("ad_numerical: Shapes of a and b must be equal.");
    }

    LA adjb = b;

    for (ptrdiff_t ii = 0; ii < p; ii++)
    {
        adjb = commutator(a, adjb);
    }

    return adjb;
}

template <typename LA>
LA ad(const LA & a, const LA & b, const int p)
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
    
    return ad_numerical<LA>(a, b, p);
}

}

#endif
