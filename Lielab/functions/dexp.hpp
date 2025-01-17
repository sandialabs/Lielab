#ifndef _LIELAB_FUNCTIONS_DEXP_HPP
#define _LIELAB_FUNCTIONS_DEXP_HPP

#include <cmath>

#include "../abstract.hpp"
#include "../domain.hpp"
#include "../utils.hpp"

#include "commutator.hpp"
#include "littlead.hpp"

namespace Lielab
{
namespace functions
{

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::gl dexp_numerical(const LA & a, const size_t order = 5)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    This is the main derivative of the exponential function. Computes it in the
    computationally intensive numerical process [1]:

    \f{equation*}{\text{dexp}_{a} = \sum_{j=0}^{\infty} \frac{1}{(j+1)!}\text{ad}_a^j \f}
    
    By default, we truncate at order 5 to best align with order 4 RK methods [2].

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    Notes
    -----

    Has a shortcut for performance.

    1. Abelian Lie algebras return identity for all orders.

    \f{equation*}{\text{dexp}_a = \sum_{j=0}^{\infty} \frac{1}{(j+1)!}\text{ad}^j_a = \text{ad}^0_a = \mathbf{I}, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    
    [2] Kenth Engø. On the construction of geometric integrators in the rkmk
        class. BIT Numerical Mathematics, 40:41–61, 2000.

    */

    Lielab::domain::gl out = ad<LA>(a, 0);

    // Special case where the domain is Abelian.
    if (a.abelian)
    {
        return out;
    }

    for (int ii = 1; ii <= order; ii++)
    {
        out = out + 1.0/factorial(ii+1)*ad<LA>(a, ii);
        // TODO: The formula on wikipedia uses -1^ii but Munthe Kaas and Engo both don't have that term?
        // I don't know how to derive this myself and not sure of any identities to check.
        // out = out + adjc*std::pow(-1.0, static_cast<double>(ii))/factorial(ii+1);
    }

    return out;
}

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::gl dexp(const LA & a, const size_t order = 5)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    Catch-all function for the derivative of the exponential function.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a Instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */

    return dexp_numerical<LA>(a, order);
}

template <Lielab::abstract::LieAlgebra LA>
LA dexp_numerical(const LA & a, const LA & b, const size_t order = 5)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    This is the main derivative of the exponential function. Computes it in the
    computationally intensive numerical process [1]:

    \f{equation*}{\text{dexp}_{a}(b) = \sum_{j=0}^{\infty} \frac{1}{(j+1)!}\text{ad}_a^j b \f}
    
    By default, we truncate at order 5 to best align with order 4 RK methods [2].

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

    \f{equation*}{\text{dexp}_a(b) = \sum_{j=0}^{\infty} \frac{1}{(j+1)!}\text{ad}^j_a(b) = \text{ad}^0_a(b) = b, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    
    [2] Kenth Engø. On the construction of geometric integrators in the rkmk
        class. BIT Numerical Mathematics, 40:41–61, 2000.

    */

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("dexp_numerical: Shapes of a and b must be equal.");
    }

    LA out(b.shape), adjc(b.shape);
    out = b;
    adjc = b;

    // Special case where the domain is abelian
    if (a.abelian)
    {
        return out;
    }

    for (int ii = 1; ii <= order; ii++)
    {
        adjc = commutator<LA>(a, adjc);
        out = out + adjc*1.0/factorial(ii+1);
        // TODO: The formula on wikipedia uses -1^ii but Munthe Kaas and Engo both don't have that term?
        // I don't know how to derive this myself and not sure of any identities to check.
        // out = out + adjc*std::pow(-1.0, static_cast<double>(ii))/factorial(ii+1);
    }

    return out;
}

template <Lielab::abstract::LieAlgebra LA>
LA dexp(const LA & a, const LA & b, const size_t order = 5)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the derivative of the exponential function.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a First instance of g
    @param[in] b Second instance of g
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of g

    */
    
    return dexp_numerical<LA>(a, b, order);
}

template <>
Lielab::domain::gl dexp(const Lielab::domain::so & x, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    Overloaded derivative of the exponential function for so.

    For \f$x \in \mathfrak{so}(3)\f$, this uses [1]:

    \f{equation*}{\theta = \Vert x \Vert, \; v = \frac{\theta}{2} \f}

    \f{equation*}{\text{dexp}(x) = \mathbf{I} + \frac{\sin^2(v)}{2v^2}\hat{x} + \frac{\theta - \sin(\theta)}{\theta^3}\hat{x}^2 \f}

    Arguments
    ---------
    @param[in] x Instance of so
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
        const Eigen::MatrixXd ad1 = Lielab::functions::ad(x, 1).get_matrix();
        const Eigen::MatrixXd ad2 = Lielab::functions::ad(x, 2).get_matrix();
        const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);

        const double theta = std::sqrt(std::pow(xbar(0), 2.0) + std::pow(xbar(1), 2.0) + std::pow(xbar(2), 2.0));
        const double v = theta/2.0;
        const double v2 = std::pow(v, 2.0);
        const double theta3 = std::pow(theta, 3.0);

        double c1 = std::pow(std::sin(v), 2.0)/(2.0*v2);
        double c2 = (theta - std::sin(theta))/(theta3);
        if (std::abs(theta) <= 1e-14)
        {
            c1 = 0.5;
            c2 = 1.0/6.0;
        }

        const Eigen::MatrixXd left = I + c1*ad1 + c2*ad2;
        return Lielab::domain::gl(left);
    }

    return dexp_numerical<Lielab::domain::so>(x, order);
}

template <>
Lielab::domain::so dexp(const Lielab::domain::so & a, const Lielab::domain::so & b, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{so} \f}
    
    Overloaded derivative of the exponential function for so.

    For \f$x \in \mathfrak{so}(2)\f$, this uses [1]:

    \f{equation*}{\text{dexp}_{a}(b) = b \f}

    For \f$x \in \mathfrak{so}(3)\f$, this uses [1]:

    \f{equation*}{\text{dexp}_{a}(b) = \text{dexp}(a)\bar{b} \f}

    For \f$a \in \mathfrak{so}(4+)\f$, this uses the numerical procedure.

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
        throw Lielab::utils::InputError("dexp: Shapes of a and b must be equal.");
    }

    if (a.shape == 2)
    {
        return b;
    }

    if (a.shape == 3)
    {
        const Lielab::domain::gl left = dexp(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::so out(3);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    return dexp_numerical<Lielab::domain::so>(a, b, order);
}

template <>
Lielab::domain::gl dexp(const Lielab::domain::se & y, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    Overloaded derivative of the exponential function for se.

    For \f$x \in \mathfrak{se}(2)\f$ (shape 3), this uses [1]:

    \f{equation*}{y = \begin{bmatrix}
    w & x \\
    \mathbf{0} & 0
    \end{bmatrix}, \; \theta = \Vert w \Vert, \; x = [u, v]^T \f}

    with

    \f{equation*}{a_\theta = \frac{\sin(\theta)}{\theta} \f}

    \f{equation*}{b_\theta = \frac{1 - \cos(\theta)}{\theta^2} \f}

    \f{equation*}{c_\theta = \frac{1 - a_\theta}{\theta} \f}

    then

    \f{equation*}{\text{dexp}(y) = \begin{bmatrix}
    a_\theta & -\theta b_\theta & c_\theta u + b_\theta v \\
    \theta b_\theta & a_\theta & c_\theta v - b_\theta u \\
    0 & 0 & 1
    \end{bmatrix} \f}

    For \f$x \in \mathfrak{se}(3)\f$ (shape 4), this uses [1]:

    \f{equation*}{y = \begin{bmatrix}
    w & x \\
    \mathbf{0} & 0
    \end{bmatrix}, \; \theta = \Vert w \Vert \f}

    with

    \f{equation*}{a_\theta = \frac{\sin(\theta)}{\theta} \f}

    \f{equation*}{b_\theta = \frac{1 - \cos(\theta)}{\theta^2} \f}

    \f{equation*}{c_\theta = \frac{1 - a_\theta}{\theta^2} \f}

    \f{equation*}{Q = \frac{a_\theta - 2 b_\theta}{\theta^2} \hat{w} + \frac{b_\theta - 3 c_\theta}{\theta^2} \hat{w}^2 \f}

    then

    \f{equation*}{\text{dexp}(y) = \begin{bmatrix}
    \text{dexp}(w) & (b_\theta \hat{x} + c_\theta (\hat{w} \hat{x} + \hat{x} \hat{w}) + (\bar{w}^T\bar{x}) Q) \\
    \mathbf{0} & \text{dexp}(w)
    \end{bmatrix} \f}

    Arguments
    ---------
    @param[in] a Instance of se.
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
        const Eigen::VectorXd ybar = y.get_vector();
        const Eigen::MatrixXd yhat = y.get_matrix();
        const Eigen::MatrixXd wmat = yhat.block(0, 0, 2, 2);
        const Eigen::MatrixXd wmat2 = wmat*wmat;
        const Eigen::VectorXd u = yhat.block(0, 2, 2, 1);

        const double theta = ybar(2);
        const double theta2 = std::pow(theta, 2.0);

        Eigen::MatrixXd left = Eigen::MatrixXd::Zero(3, 3);

        // SO component.
        double atheta = std::sin(theta)/theta;
        double btheta = (1.0 - std::cos(theta))/theta2;
        // Note: I think Eade has a bug in Eq. 83. He used ctheta = (1.0 - atheta)/(theta^2)
        //       but I think it should be ctheta = (1.0 - atheta)/theta. This is in his
        //       November 18, 2018 version of the document. I'm calling this dtheta for clarity.
        double dtheta = (1.0 - atheta)/theta;
        if (std::abs(theta) <= 1e-14)
        {
            atheta = 1.0;
            btheta = 0.5;
            dtheta = 0.0;
        }
        left(0,0) = atheta;
        left(0,1) = -theta*btheta;
        left(1,0) = theta*btheta;
        left(1,1) = atheta;

        // R component.
        left(0,2) = dtheta*u(0) + btheta*u(1);
        left(1,2) = dtheta*u(1) - btheta*u(0);
        left(2,2) = 1.0;
        return Lielab::domain::gl(left);
    }

    if (y.shape == 4)
    {
        // Source: Eade
        const Eigen::VectorXd ybar = y.get_vector();
        const Eigen::MatrixXd yhat = y.get_matrix();
        const Eigen::MatrixXd wmat = yhat.block(0, 0, 3, 3);
        const Eigen::MatrixXd wmat2 = wmat*wmat;
        const Eigen::VectorXd u = yhat.block(0, 3, 3, 1);
        Lielab::domain::so ux(3);
        ux.set_vector(u);
        const Eigen::MatrixXd umat = ux.get_matrix();

        const double wmag2 = std::pow(ybar(3), 2.0) + std::pow(ybar(4), 2.0) + std::pow(ybar(5), 2.0);
        const double wmag = std::sqrt(wmag2);

        Eigen::MatrixXd left = Eigen::MatrixXd::Zero(6, 6);

        // SO component. Re-use SO calculation.
        const Lielab::domain::so w(wmat);
        const Lielab::domain::gl dw = Lielab::functions::dexp(w);
        left.block(0, 0, 3, 3) = dw.get_matrix();
        left.block(3, 3, 3, 3) = dw.get_matrix();

        // R component.
        const Eigen::VectorXd wvec = w.get_vector();
        double aw = std::sin(wmag)/wmag;
        double bw = (1.0 - std::cos(wmag))/wmag2;
        double cw = (1.0 - aw)/wmag2;
        double q1 = (aw - 2.0*bw)/wmag2;
        double q2 = (bw - 3.0*cw)/wmag2;
        if (std::abs(wmag) <= 1e-14)
        {
            aw = 1.0;
            bw = 0.5;
            cw = 1.0/6.0;
            q1 = -1.0/12.0;
            q2 = -1.0/60.0;
        }
        const Eigen::MatrixXd Q = q1*wmat + q2*wmat2;
        left.block(0, 3, 3, 3) = (bw*umat + cw*(wmat*umat + umat*wmat) + (wvec(0)*u(0) + wvec(1)*u(1) + wvec(2)*u(2))*Q);
        return Lielab::domain::gl(left);
    }

    return dexp_numerical<Lielab::domain::se>(y, order);
}

template <>
Lielab::domain::se dexp(const Lielab::domain::se & a, const Lielab::domain::se & b, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
    
    Overloaded derivative of the exponential function for se.

    For \f$a \in \mathfrak{se}(2)\f$ (shape 3), this uses [1]:

    \f{equation*}{\text{dexp}_{a}(b) = \text{dexp}(a)\bar{b} \f}

    For \f$a \in \mathfrak{se}(3)\f$ (shape 4), this uses [1]:

    \f{equation*}{\text{dexp}_{a}(b) = \text{dexp}(a)\bar{b} \f}

    For \f$a \in \mathfrak{se}(4+)\f$ (shape 5+), this uses the numerical procedure.

    Arguments
    ---------
    @param[in] a First instance of se.
    @param[in] b Second instance of se.
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of se.

    References
    ----------
    [1] Ethan Eade. Derivative of the exponential map. Technical report,
        November 2018. [Online]. Available: https://ethaneade.com/exp_diff.pdf.

    */

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("dexp: Shapes of a and b must be equal.");
    }

    if (a.shape == 3)
    {
        const Lielab::domain::gl left = dexp(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::se out(3);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    if (a.shape == 4)
    {
        const Lielab::domain::gl left = dexp(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::se out(4);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    return dexp_numerical<Lielab::domain::se>(a, b, order);
}

template <>
Lielab::domain::gl dexp(const Lielab::domain::su & x, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    Overloaded derivative of the exponential function for su.

    For \f$x \in \mathfrak{su}(2)\f$, this uses:

    // \f{equation*}{\theta = \Vert x \Vert, \; v = \frac{\theta}{2} \f}

    // \f{equation*}{\text{dexp}(x) = \mathbf{I} + \frac{\sin^2(v)}{2v^2}\hat{x} + \frac{\theta - \sin(\theta)}{\theta^3}\hat{x}^2 \f}

    Arguments
    ---------
    @param[in] x Instance of su
    @param[in] order Order of the series expansion.
    @param[out] out An instance of gl

    References
    ----------
    Derived it myself - Mike Sparapany

    */

    // if (x.shape == 1)
    // {
    //     return Lielab::domain::gl(Eigen::MatrixXd::Identity(1,1));
    // }

    if (x.shape == 2)
    {
        const Eigen::Vector3d xbar = x.get_vector();
        const Eigen::MatrixXd ad1 = Lielab::functions::ad(x, 1).get_matrix();
        const Eigen::MatrixXd ad2 = Lielab::functions::ad(x, 2).get_matrix();
        const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);

        const double theta = (std::sqrt(std::pow(xbar(0), 2.0) + std::pow(xbar(1), 2.0) + std::pow(xbar(2), 2.0)));
        const double theta2 = std::pow(theta, 2.0);
        const double theta3 = std::pow(theta, 3.0);

        double c1 = std::pow(std::sin(theta), 2.0)/(2.0*theta2);
        double c2 = (2.0*theta - std::sin(2.0*theta))/(8.0*theta3);
        if (std::abs(theta) <= 1e-14)
        {
            c1 = 0.5;
            c2 = 1.0/6.0;
        }

        const Eigen::MatrixXd left = I + c1*ad1 + c2*ad2;
        return Lielab::domain::gl(left);
    }

    return dexp_numerical<Lielab::domain::su>(x, order);
}

template <>
Lielab::domain::su dexp(const Lielab::domain::su & a, const Lielab::domain::su & b, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    
    Overloaded derivative of the exponential function for su.

    For \f$x \in \mathfrak{su}(2)\f$, this uses:

    \f{equation*}{\text{dexp}_{a}(b) = \text{dexp}(a)\bar{b} \f}

    For \f$a \in \mathfrak{su}(3+)\f$ and \f$\mathfrak{su}(2)\f$, this uses the numerical procedure.

    Arguments
    ---------
    @param[in] a First instance of so
    @param[in] b Second instance of so
    @param[in] order Order of the series expansion.
    @param[out] out An instance of so

    References
    ----------
    Derived it myself - Mike Sparapany

    */

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("dexp: Shapes of a and b must be equal.");
    }

    // if (a.shape == 1)
    // {
    //     return b;
    // }

    if (a.shape == 2)
    {
        const Lielab::domain::gl left = dexp(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::su out(2);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    return dexp_numerical<Lielab::domain::su>(a, b, order);
}

}
}

#endif
