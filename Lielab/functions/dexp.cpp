#include "dexp.hpp"

#include "commutator.hpp"
#include "littlead.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/utils.hpp"

#include <cmath>

namespace Lielab::functions
{

template <>
Lielab::domain::glr dexp(const Lielab::domain::glr & x, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overloaded derivative of the exponential function for glr.

    Arguments
    ---------
    @param[in] x Instance of glr
    @param[in] order Order of the series expansion.
    @param[out] out An instance of glr

    References
    ----------
    [1] I made it up.

    */

    const size_t shape = x.get_shape();

    if (shape == 2)
    {
        // Ref [1]

        const Eigen::Matrix2d xhat = x.get_matrix();
        const double a = xhat(0,0);
        const double b = xhat(0,1);
        const double c = xhat(1,0);
        const double d = xhat(1,1);
        const double theta = std::pow(a - d, 2.0) + 4*b*c;

        const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(4, 4);
        const Eigen::MatrixXd ad1 = Lielab::functions::ad(x, 1).get_matrix();
        const Eigen::MatrixXd ad2 = Lielab::functions::ad(x, 2).get_matrix();
        
        double c1;
        double c2;

        if (std::abs(theta) <= 1e-14)
        {
            c1 = 0.5;
            c2 = 1.0/6.0;
        }
        else if (theta <= -1e-14)
        {
            c1 = (std::cos(std::sqrt(-theta)) - 1.0)/theta;
            c2 = -(std::sin(std::sqrt(-theta)) - std::sqrt(-theta))/std::pow(-theta, 3.0/2.0);
        }
        else
        {
            c1 = (std::cosh(std::sqrt(theta)) - 1.0)/theta;
            c2 = (std::sinh(std::sqrt(theta)) - std::sqrt(theta))/std::pow(theta, 3.0/2.0);
        }

        return Lielab::domain::glr(I + c1*ad1 + c2*ad2);
    }

    return dexp_numerical<Lielab::domain::glr>(x, order);
}

template <>
Lielab::domain::glr dexp(const Lielab::domain::glr & x, const Lielab::domain::glr & y, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overloaded derivative of the exponential function for glr.

    Arguments
    ---------
    @param[in] a First instance of glr
    @param[in] b Second instance of glr
    @param[in] order Order of the series expansion.
    @param[out] out An instance of glr

    References
    ----------
    [1] Derived it myself.

    */
    
    const size_t shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dexp: Shapes of x and y must be equal.");
    }

    if (shape == 2)
    {
        // Ref [1]

        const Eigen::MatrixXd adxhat = dexp<Lielab::domain::glr>(x).get_matrix();
        const Eigen::MatrixXd ybar = y.get_vector();
        Lielab::domain::glr out(2);
        out.set_vector(adxhat*ybar);
        return out;
    }

    return dexp_numerical<Lielab::domain::glr>(x, y, order);
}

template <>
Lielab::domain::glr dexp(const Lielab::domain::so & x, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overloaded derivative of the exponential function for so.

    For \f$x \in \mathfrak{so}(3)\f$, this uses [1]:

    \f{equation*}{\theta = \Vert x \Vert, \; v = \frac{\theta}{2} \f}

    \f{equation*}{\text{dexp}(x) = \mathbf{I} + \frac{\sin^2(v)}{2v^2}\hat{x} + \frac{\theta - \sin(\theta)}{\theta^3}\hat{x}^2 \f}

    Arguments
    ---------
    @param[in] x Instance of so
    @param[in] order Order of the series expansion.
    @param[out] out An instance of glr

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.

    */

    const size_t shape = x.get_shape();

    if (shape == 2)
    {
        return Lielab::domain::glr(Eigen::MatrixXd::Identity(1,1));
    }

    if (shape == 3)
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
        return Lielab::domain::glr(left);
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

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("dexp: Shapes of a and b must be equal.");
    }

    if (shape == 2)
    {
        return b;
    }

    if (shape == 3)
    {
        const Lielab::domain::glr left = dexp(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::so out(3);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    return dexp_numerical<Lielab::domain::so>(a, b, order);
}

template <>
Lielab::domain::glr dexp(const Lielab::domain::se & y, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
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
    @param[out] out An instance of glr

    References
    ----------
    [1] Ethan Eade. Derivative of the exponential map. Technical report,
        November 2018. [Online]. Available: https://ethaneade.com/exp_diff.pdf.

    */

    const size_t shape = y.get_shape();

    if (shape == 3)
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
        return Lielab::domain::glr(left);
    }

    if (shape == 4)
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
        const Lielab::domain::glr dw = Lielab::functions::dexp(w);
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
        return Lielab::domain::glr(left);
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

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("dexp: Shapes of a and b must be equal.");
    }

    if (shape == 3)
    {
        const Lielab::domain::glr left = dexp(a);
        const Eigen::MatrixXd y = b.get_vector();
        return Lielab::domain::se::from_vector(left.get_matrix()*y);
    }

    if (shape == 4)
    {
        const Lielab::domain::glr left = dexp(a);
        const Eigen::MatrixXd y = b.get_vector();
        return Lielab::domain::se::from_vector(left.get_matrix()*y);
    }

    return dexp_numerical<Lielab::domain::se>(a, b, order);
}

template <>
Lielab::domain::glr dexp(const Lielab::domain::su & x, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overloaded derivative of the exponential function for su.

    For \f$x \in \mathfrak{su}(2)\f$, this uses:

    // \f{equation*}{\theta = \Vert x \Vert, \; v = \frac{\theta}{2} \f}

    // \f{equation*}{\text{dexp}(x) = \mathbf{I} + \frac{\sin^2(v)}{2v^2}\hat{x} + \frac{\theta - \sin(\theta)}{\theta^3}\hat{x}^2 \f}

    Arguments
    ---------
    @param[in] x Instance of su
    @param[in] order Order of the series expansion.
    @param[out] out An instance of glr

    References
    ----------
    Derived it myself - Mike Sparapany

    */

    const size_t shape = x.get_shape();

    // if (shape == 1)
    // {
    //     return Lielab::domain::glr(Eigen::MatrixXd::Identity(1,1));
    // }

    if (shape == 2)
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
        return Lielab::domain::glr(left);
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

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("dexp: Shapes of a and b must be equal.");
    }

    // if (shape == 1)
    // {
    //     return b;
    // }

    if (shape == 2)
    {
        const Lielab::domain::glr left = dexp(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::su out(2);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    return dexp_numerical<Lielab::domain::su>(a, b, order);
}

template <>
Lielab::domain::CompositeAlgebra dexp(const Lielab::domain::CompositeAlgebra & a, const Lielab::domain::CompositeAlgebra & b, const size_t order)
{
    /*!
    * CompositeAlgebra derivative of the exponential overload.
    */

    Lielab::domain::CompositeAlgebra out;

    for (size_t ii = 0; ii < a.space.size(); ii++)
    {
        const size_t ind = a.space[ii].index();
        if (ind == Lielab::domain::CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::cn>(a.space[ii]),
                                                        std::get<Lielab::domain::cn>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::glr>(a.space[ii]),
                                                        std::get<Lielab::domain::glr>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::glc>(a.space[ii]),
                                                        std::get<Lielab::domain::glc>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::rn>(a.space[ii]),
                                                        std::get<Lielab::domain::rn>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::se>(a.space[ii]),
                                                        std::get<Lielab::domain::se>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::so>(a.space[ii]),
                                                        std::get<Lielab::domain::so>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::sp>(a.space[ii]),
                                                        std::get<Lielab::domain::sp>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::su>(a.space[ii]),
                                                        std::get<Lielab::domain::su>(b.space[ii]), order));
        }
    }

    return out;
}

}
