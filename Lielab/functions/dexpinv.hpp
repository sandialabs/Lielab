#ifndef _LIELAB_FUNCTIONS_DEXPINV_HPP
#define _LIELAB_FUNCTIONS_DEXPINV_HPP

#include "../abstract.hpp"
#include "../domain.hpp"

#include "dexp.hpp"

namespace Lielab
{
namespace functions
{

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::gl dexpinv_numerical(const LA & a, const size_t order = 5)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    This is the main inverse derivative of the exponential function. Computes
    it in the computationally intensive numerical process [1]:

    \f{equation*}{\text{dexp}_{a}^{-1} = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}^j_a \f}
    
    Where \f$B_j\f$ are Bernoulli numbers. By default, we truncate at order 5
    to best align with order 4 RK methods [2].

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    Notes
    -----

    Has a shortcut for performance.

    1. Abelian Lie algebras return identity for all orders.

    \f{equation*}{\text{dexpinv}_a = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}^j_a = \text{ad}^0_a = \mathbf{I}, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    
    [2] Kenth Engø. On the construction of geometric integrators in the rkmk
        class. BIT Numerical Mathematics, 40:41–61, 2000.

    */

    Lielab::domain::gl out = ad<LA>(a, 0);

    // Special case where the domain is abelian
    if (a.abelian)
    {
        return out;
    }

    for (int ii = 1; ii <= order; ii++)
    {
        out = out + bernoulli(ii)/factorial(ii)*ad<LA>(a, ii);
    }

    return out;
}

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::gl dexpinv(const LA & a, const size_t order = 5)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the inverse derivative of the exponential function.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a Instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexpinv_numerical<LA>(a, order);
}

template <Lielab::abstract::LieAlgebra LA>
LA dexpinv_numerical(const LA & a, const LA & b, const size_t order = 5)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    This is the main inverse derivative of the exponential function. Computes
    it in the computationally intensive numerical process [1]:

    \f{equation*}{\text{dexp}_{a}^{-1}(b) = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}_a^j b \f}
    
    Where \f$B_j\f$ are Bernoulli numbers. By default, we truncate at order 5
    to best align with order 4 RK methods [2].

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] b Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    Notes
    -----

    Has a shortcut for performance.

    1. Abelian Lie algebras return b for all orders.

    \f{equation*}{\text{dexpinv}_a(b) = \sum_{j=0}^{\infty} \frac{B_j}{j!}\text{ad}^j_a(b) = \text{ad}^0_a(b) = b, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    
    [2] Kenth Engø. On the construction of geometric integrators in the rkmk
        class. BIT Numerical Mathematics, 40:41–61, 2000.

    */

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("dexpinv_numerical: Shapes of a and b must be equal.");
    }
    
    LA out(b.shape), adjc(b.shape);
    adjc = b;
    out = adjc;

    // Special case where the domain is abelian
    if (a.abelian)
    {
        return out;
    }

    for (int ii = 1; ii <= order; ii++)
    {
        adjc = commutator<LA>(a, adjc);
        out = out + adjc*bernoulli(ii)/factorial(ii);
    }

    return out;
}

template <Lielab::abstract::LieAlgebra LA>
LA dexpinv(const LA & a, const LA & b, const size_t order = 5)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the inverse derivative of the exponential function.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] b Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexpinv_numerical<LA>(a, b, order);
}

template <>
Lielab::domain::gl dexpinv(const Lielab::domain::so & x, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    Overloaded inverse derivative of the exponential function for so.

    For \f$x \in \mathfrak{so}(3)\f$, this uses [1]:

    \f{equation*}{\theta = \Vert x \Vert, \; v = \frac{\theta}{2} \f}

    \f{equation*}{\text{dexp}^{-1}(x) = \mathbf{I} - \frac{1}{2}\hat{x} - \frac{\theta \cot(v) - 2}{2 \theta^2}\hat{x}^2 \f}

    Arguments
    ---------
    @param[in] x First instance of so
    @param[in] order Order of the series expansion.
    @param[out] out An instance of gl

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.

    */

    if (x.shape == 2)
    {
        return Lielab::domain::gl(Eigen::MatrixXd::Identity(1,1));
    }

    if (x.shape == 3)
    {
        // Source: Iserles
        const Eigen::Vector3d xbar = x.get_vector();
        const Eigen::MatrixXd xhat = x.get_matrix();
        const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);

        const double theta = std::sqrt(std::pow(xbar(0), 2.0) + std::pow(xbar(1), 2.0) + std::pow(xbar(2), 2.0));
        const double v = theta/2.0;
        const double sv = std::sin(v);
        const double cv = std::cos(v);
        const double theta2 = std::pow(theta, 2.0);
        const Eigen::MatrixXd xhat2 = xhat*xhat;

        double c2 = (theta*cv/sv - 2.0)/(2*theta2);
        if (std::abs(theta) <= 1e-14)
        {
            c2 = -1.0/12.0;
        }

        const Eigen::MatrixXd left = I - 1.0/2.0*xhat - c2*xhat2;
        return Lielab::domain::gl(left);
    }

    return dexpinv_numerical<Lielab::domain::so>(x, order);
}

template <>
Lielab::domain::so dexpinv(const Lielab::domain::so & a, const Lielab::domain::so & b, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{so} \f}
    
    Overloaded inverse derivative of the exponential function for so.

    For \f$x \in \mathfrak{so}(3)\f$, this uses [1]:

    \f{equation*}{\text{dexp}^{-1}_{x}(y) = \text{dexp}^{-1}(x)\bar{y} \f}

    For \f$a \in \mathfrak{so}(4+)\f$ and \f$\mathfrak{so}(2)\f$, this uses the numerical procedure.

    Arguments
    ---------
    @param[in] a First instance of so
    @param[in] b Second instance of so
    @param[in] order Order of the series expansion.
    @param[out] out An instance of so

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.

    */

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("dexpinv: Shapes of a and b must be equal.");
    }

    if (a.shape == 2)
    {
        return b;
    }

    if (a.shape == 3)
    {
        const Lielab::domain::gl left = dexpinv(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::so out(3);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    return dexpinv_numerical<Lielab::domain::so>(a, b, order);
}

template <>
Lielab::domain::gl dexpinv(const Lielab::domain::se & y, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    Overloaded inverse derivative of the exponential function for se.

    For \f$y \in \mathfrak{se}(2)\f$ (shape 3), this uses [1]:

    \f{equation*}{\text{dexp}(y) = \begin{bmatrix}
    A & b \\
    \mathbf{0} & 1
    \end{bmatrix} \f}

    then

    \f{equation*}{\text{dexp}^{-1}(y) = \begin{bmatrix}
    A^{-1} & -A^{-1} b \\
    \mathbf{0} & 1
    \end{bmatrix} \f}

    For \f$y \in \mathfrak{se}(3)\f$ (shape 4), this uses [1]:

    \f{equation*}{\text{dexp}(y) = \begin{bmatrix}
    A & B \\
    \mathbf{0} & A
    \end{bmatrix} \f}

    then

    \f{equation*}{\text{dexp}^{-1}(y) = \begin{bmatrix}
    A^{-1} & -A^{-1} B A^{-1} \\
    \mathbf{0} & A^{-1}
    \end{bmatrix} \f}

    Arguments
    ---------
    @param[in] a Instance of se
    @param[in] order Order of the series expansion.
    @param[out] out An instance of gl

    References
    ----------
    [1] Ethan Eade. Derivative of the exponential map. Technical report,
        November 2018. [Online]. Available: https://ethaneade.com/exp_diff.pdf.

    */

    if (y.shape == 3)
    {
        // Source: Eade
        const Lielab::domain::gl left_inv = Lielab::functions::dexp(y);
        const Eigen::MatrixXd A = left_inv.get_matrix().block(0, 0, 2, 2);
        const Eigen::VectorXd v = left_inv.get_matrix().block(0, 2, 2, 1);
        const Eigen::MatrixXd Ainv = A.inverse();

        Eigen::MatrixXd left = Eigen::MatrixXd::Zero(3, 3);
        left.block(0, 0, 2, 2) = Ainv;
        left.block(0, 2, 2, 1) = -Ainv*v;
        left(2, 2) = 1.0;
        return Lielab::domain::gl(left);
    }

    if (y.shape == 4)
    {
        // Source: Eade
        const Lielab::domain::gl left_inv = Lielab::functions::dexp(y);
        const Eigen::MatrixXd A = left_inv.get_matrix().block(0, 0, 3, 3);
        const Eigen::MatrixXd B = left_inv.get_matrix().block(0, 3, 3, 3);
        const Eigen::MatrixXd Ainv = A.inverse();

        Eigen::MatrixXd left = Eigen::MatrixXd::Zero(6, 6);

        left.block(0, 0, 3, 3) = Ainv;
        left.block(0, 3, 3, 3) = -Ainv*B*Ainv;
        left.block(3, 3, 3, 3) = Ainv;
        return Lielab::domain::gl(left);
    }

    return dexpinv_numerical<Lielab::domain::se>(y, order);
}

template <>
Lielab::domain::se dexpinv(const Lielab::domain::se & a, const Lielab::domain::se & b, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
    
    Overloaded derivative of the exponential function for se.

    For \f$a \in \mathfrak{se}(2)\f$ (shape 3), this uses [1]:

    \f{equation*}{\text{dexp}^{-1}_{a}(b) = \text{dexp}^{-1}(a)\bar{b} \f}

    For \f$a \in \mathfrak{se}(3)\f$ (shape 4), this uses [1]:

    \f{equation*}{\text{dexp}^{-1}_{a}(b) = \text{dexp}^{-1}(a)\bar{b} \f}

    For \f$a \in \mathfrak{se}(4+)\f$ (shape 5+), this uses the numerical procedure.

    Arguments
    ---------
    @param[in] a First instance of se
    @param[in] b Second instance of se
    @param[in] order Order of the series expansion.
    @param[out] out An instance of se

    References
    ----------
    [1] Ethan Eade. Derivative of the exponential map. Technical report,
        November 2018. [Online]. Available: https://ethaneade.com/exp_diff.pdf.

    */

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("dexpinv: Shapes of a and b must be equal.");
    }

    if (a.shape == 3)
    {
        const Lielab::domain::gl left = dexpinv(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::se out(3);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    if (a.shape == 4)
    {
        const Lielab::domain::gl left = dexpinv(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::se out(4);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    return dexp_numerical<Lielab::domain::se>(a, b, order);
}

}
}

#endif
