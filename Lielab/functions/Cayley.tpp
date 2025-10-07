#ifndef LIELAB_FUNCTIONS_CAYLEY_TPP
#define LIELAB_FUNCTIONS_CAYLEY_TPP

#include "Cayley.hpp"

#include "commutator.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <typename LA>
Lielab::domain::LieIII<LA> Cayley(const LA & g)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
    *
    * Cayley transform.
    * 
    * \f{equation*}{ Cayley(g) = \frac{Id_g + g/2}{Id_g - g/2} \f}
    * 
    * @param[in] g A Lie algebra.
    * @param[out] G A Lie group.
    * 
    * Source: Engø, Kenth. "On the construction of geometric integrators in the RKMK class."
    * BIT Numerical Mathematics 40.1 (2000): 41-61.
    * 
    * TODO: Restrict this to only work with O, SO, and SP.
    * TODO: Check that this works with SU. Math says no but simulations say otherwise.
    */

    const size_t shape = g.get_shape();

    const Eigen::MatrixXd m = g.get_matrix();
    const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(shape, shape);

    return (Id + m/2.0)*(Id - m/2.0).inverse();
}

template <typename LA>
Lielab::domain::LieIII<LA> Cayley2(const LA & g)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
    *
    * Cayley transform.
    * 
    * \f{equation*}{ Cayley2(g) = Cayley(\xi_1 x_1)Cayley(\xi_2 x_2) \cdots Cayley(\xi_n x_n) \f}
    * 
    * @param[in] g A Lie algebra.
    * @param[out] G A Lie group.
    * 
    * Source: Engø, Kenth. "On the construction of geometric integrators in the RKMK class."
    * BIT Numerical Mathematics 40.1 (2000): 41-61.
    * 
    * TODO: Check that this works with SU. Math says no but simulations say otherwise.
    */

    const size_t dim = g.get_dimension();
    const size_t shape = g.get_shape();
    const Eigen::VectorXd v = g.get_vector();

    Lielab::domain::LieIII<LA> out(shape);

    for (size_t ii = 0; ii < dim; ii++)
    {
        const LA h = LA::basis(ii, shape);
        out = out*Cayley(v[ii]*h);
    }

    return out;
}

template <typename LA>
LA dCayleyinv(const LA & u, const LA & v)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}) \rightarrow \mathfrak{g} \f}
    * 
    * Derivative of the inverse Cayley function.
    *
    * \f{equation*}{ \text{dCayley}_u^{-1}(v) = v - \frac{1}{2} [u,v] - \frac{1}{4} u \cdot v \cdot u \f}
    * 
    * Source: Engø, Kenth. "On the construction of geometric integrators in the RKMK class."
    * BIT Numerical Mathematics 40.1 (2000): 41-61.
    */

    const size_t shape = u.get_shape();

    if (shape != v.get_shape())
    {
        throw Lielab::utils::InputError("dCayleyinv: Shapes of a and b must be equal.");
    }

    const typename LA::matrix_t uhat = u.get_matrix();
    const typename LA::matrix_t vhat = v.get_matrix();

    LA temp(shape); // TODO: Make this ::from_shape
    temp = v;

    if (!v.abelian)
    {
        temp = temp - 1.0/2.0*commutator(u, v);
    }

    return temp - LA(1.0/4.0*uhat*vhat*uhat);
}

}

#endif
