#ifndef LIELAB_FUNCTIONS_CAYLEY_TPP
#define LIELAB_FUNCTIONS_CAYLEY_TPP

#include "Cayley.hpp"

#include "adjoint.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <typename g>
Lielab::domain::LieIII<g> cay(const g& x)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
    *
    * Cayley transform.
    * 
    * \f{equation*}{ cay(x) = \frac{Id + \hat{x}/2}{Id - \hat{x}/2} \f}
    * 
    * This is only valid for skew Hermitian algebras to guarantee the inverse
    * exists, but this is not checked and therefore this function can be,
    * potentially improperly, called with any Lie algebra.
    * 
    * @param[in] x A Lie algebra.
    * @param[out] g A Lie group.
    * 
    * References
    * ----------
    *     [1] Engø, Kenth. "On the construction of geometric integrators in the RKMK class."
    *         BIT Numerical Mathematics 40.1 (2000): 41-61.
    */

    const size_t shape = x.get_shape();

    const typename g::matrix_t xhat = x.get_matrix();
    const typename g::matrix_t Id = g::matrix_t::Identity(shape, shape);

    // TODO: Abelian speedhack, just add Identity and return??

    return typename g::matrix_t((Id + xhat/2.0)*(Id - xhat/2.0).inverse());
}

template <typename G>
Lielab::domain::LieIII<G> cayinv(const G& g)
{
    /*! \f{equation*}{ (G) \rightarrow \mathfrak{g} \f}
    *
    * Inverse Cayley transform. Non-inverse formula from Ref 1, then
    * I just inverted it.
    * 
    * This is only valid for unitary groups to guarantee the inverse
    * exists, but this is not checked and therefore this function can be,
    * potentially improperly, called with any Lie group.
    * 
    * @param[in] g A Lie group.
    * @param[out] x A Lie algebra.
    * 
    * References
    * ----------
    *     [1] Engø, Kenth. "On the construction of geometric integrators in the RKMK class."
    *         BIT Numerical Mathematics 40.1 (2000): 41-61.
    *     [2] Inverted the equation myself.
    */

    using Lielab::domain::LieIII;

    const size_t shape = g.get_shape();

    const typename G::matrix_t ghat = g.get_matrix();
    const typename G::matrix_t Id = G::matrix_t::Identity(shape, shape);

    // TODO: Abelian speedhack, just subtract Identity and return??

    return typename LieIII<G>::matrix_t(2.0*((ghat + Id).inverse()*(ghat - Id)));
}

// template <typename LA>
// Lielab::domain::glr dcay(const LA& x)
// {
// This is wrong   
//     const size_t shape = x.get_shape();
//     const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(shape, shape);
//     const Eigen::MatrixXd adxhat = ad(x).get_matrix();

//     return (Id - 0.5*adxhat).inverse()*(Id + 0.5*adxhat).inverse();
//     // return ((2.0*Id - adxhat).pow(-1).matrix() + (2.0*Id + adxhat)*(2.0*Id - adxhat).pow(-2).matrix());
// }

template <typename g>
g dcay(const g& x, const g& y)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    *
    * Differential of the Cayley transform. Ref 1 Equation 3.5 scaled by 2.
    * 
    * @param[in] x A Lie algebra.
    * @param[in] y A Lie algebra.
    * @param[out] z A Lie algebra.
    * 
    * References
    * ----------
    *     [1] Müller, Andreas. "Review of the exponential and Cayley map on SE (3) as relevant for Lie group
    *         integration of the generalized Poisson equation and flexible multibody systems." Proceedings of
    *         The Royal Society A 477.2253 (2021): 20210303.
    */

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dcay: Shapes of x and y must be equal.");
    }

    const typename g::matrix_t xhat = x.get_matrix();
    const typename g::matrix_t yhat = y.get_matrix();
    const typename g::matrix_t Id = g::matrix_t::Identity(shape, shape);

    return typename g::matrix_t((Id - xhat/2.0).inverse()*yhat*(Id + xhat/2.0).inverse());
}

// template <typename LA>
// Lielab::domain::glr dcayinv(const LA& x)
// {
// This is also wrong
//     /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
//     *
//     * Differential inverse of the Cayley transform.
//     * 
//     * @param[in] g A Lie algebra.
//     * @param[out] G A Lie group.
//     * 
//     * References
//     * ----------
//     *     [1] Engø, Kenth. "On the construction of geometric integrators in the RKMK class."
//     *         BIT Numerical Mathematics 40.1 (2000): 41-61.
//     */
    
//     const size_t shape = x.get_shape();
//     const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(shape, shape);
//     const Eigen::MatrixXd adxhat = ad(x).get_matrix();
//     const Eigen::VectorXd xbar = x.get_vector();

//     // return (Id + 0.5*adxhat)*(Id - 0.5*adxhat);
//     return (Id - 0.5*adxhat + 0.25*xbar*xbar.transpose());
// }

template <typename g>
g dcayinv(const g& x, const g& y)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    
    Derivative of the inverse Cayley function.
    
    \f{equation*}{ \text{dcay}_x^{-1}(y) = y - \frac{1}{2} [x,y] - \frac{1}{4} x \cdot y \cdot x \f}
    
    References
    ----------
    Engø, Kenth. "On the construction of geometric integrators in the RKMK class."
    BIT Numerical Mathematics 40.1 (2000): 41-61.
    */

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dcayinv: Shapes of x and y must be equal.");
    }

    const typename g::matrix_t xhat = x.get_matrix();
    const typename g::matrix_t yhat = y.get_matrix();

    typename g::matrix_t temp = yhat;

    if (!y.is_abelian())
    {
        temp -= 1.0/2.0*(xhat*yhat - yhat*xhat);
    }

    return typename g::matrix_t(temp - 1.0/4.0*xhat*yhat*xhat);
}

template <typename g>
Lielab::domain::LieIII<g> cay2(const g& x)
{
    /*! \f{equation*}{ (\mathfrak{x}) \rightarrow G \f}
    *
    * Cayley transform.
    * 
    * \f{equation*}{ cay2(x) = cay(\xi_1 x_1) cay(\xi_2 x_2) \cdots cay(\xi_n x_n) \f}
    * 
    * @param[in] x A Lie algebra.
    * @param[out] g A Lie group.
    * 
    * Source: Engø, Kenth. "On the construction of geometric integrators in the RKMK class."
    * BIT Numerical Mathematics 40.1 (2000): 41-61.
    */

    const int dim = x.get_dimension();
    const int shape = x.get_shape();
    const Eigen::VectorXd xbar = x.get_vector();

    Lielab::domain::LieIII<g> out(shape); // TODO: use ::zero instead.

    for (int ii = 0; ii < dim; ii++)
    {
        const g h = g::basis(ii, shape);
        out *= cay(xbar(ii)*h);
    }

    return out;
}

}

#endif
