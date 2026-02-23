#include "exponential.hpp"

#include "adjoint.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/testing/assertions.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <format>
#include <tuple>

namespace Lielab::functions
{

template <>
Lielab::domain::CN exp(const Lielab::domain::cn& x)
{
    /*! \f{equation*}{ (\mathfrak{cn}) \rightarrow CN \f}
    
    Exponential function overload for \f$\mathfrak{cn}\f$.

    \f$\mathfrak{cn}\f$ is Abelian so the data is copied directly.
    
    Arguments
    ---------
    @param[in] x An instance of cn
    @param[out] out An instance of CN

    */

    return Lielab::domain::CN::from_complex_vector(x.point);
}

template <>
Lielab::domain::GLC exp(const Lielab::domain::glc& x)
{
    /*! \f{equation*}{ (\mathfrak{glc}) \rightarrow GLC \f}
    
    Exponential function overload for \f$\mathfrak{glc}\f$.

    Arguments
    ---------
    @param[in] x An instance of glc
    @param[out] out An instance of GLC

    References
    ----------
    [1] Bernstein, Dennis S., and Wasin So. "Some explicit formulas for the matrix
        exponential." IEEE Transactions on Automatic Control 38.8 (1993): 1228-1232.

    */

    const size_t shape = x.get_shape();

    if (shape == 2)
    {
        // Ref [1]
        const double eps = std::numeric_limits<double>::epsilon();

        const Eigen::Matrix2cd xhat = x.get_matrix();
        Eigen::Matrix2cd Xhat = Eigen::Matrix2cd::Zero(2,2);

        const std::complex<double> a = xhat(0,0);
        const std::complex<double> b = xhat(0,1);
        const std::complex<double> c = xhat(1,0);
        const std::complex<double> d = xhat(1,1);

        const std::complex<double> delta2 = std::pow(a - d, 2.0) + std::complex<double>(4.0, 0.0)*b*c;
        const std::complex<double> m = std::exp((a + d)/2.0);

        if (std::abs(delta2) < 2.0*eps)
        {
            // Ref [1] Corollary 2.3 Case 2.
            Xhat(0,0) = m*(1.0 + (a - d)/2.0);
            Xhat(0,1) = m*(b);
            Xhat(1,0) = m*(c);
            Xhat(1,1) = m*(1.0 - (a - d)/2.0);
        }
        else
        {
            // Ref [1] Corollary 2.3 Case 1.
            const std::complex<double> delta = std::complex<double>(0.5, 0.0)*std::sqrt(delta2);
            Xhat(0,0) = m*(std::cosh(delta) + (a - d)/2.0*std::sinh(delta)/delta);
            Xhat(0,1) = m*(b*std::sinh(delta)/delta);
            Xhat(1,0) = m*(c*std::sinh(delta)/delta);
            Xhat(1,1) = m*(std::cosh(delta) - (a - d)/2.0*std::sinh(delta)/delta);
        }
        
        return Lielab::domain::GLC(Xhat);
    }

    return exp_numerical(x);
}

template <>
Lielab::domain::GLR exp(const Lielab::domain::glr& x)
{
    /*! \f{equation*}{ (\mathfrak{glr}) \rightarrow GLR \f}
    
    Exponential function overload for \f$\mathfrak{glr}\f$.

    Arguments
    ---------
    @param[in] x An instance of glr
    @param[out] out An instance of GLR

    References
    ----------
    [1] Bernstein, Dennis S., and Wasin So. "Some explicit formulas for the matrix
        exponential." IEEE Transactions on Automatic Control 38.8 (1993): 1228-1232.

    */

    const size_t shape = x.get_shape();

    if (shape == 2)
    {
        // Ref [1]
        const double eps = std::numeric_limits<double>::epsilon();

        const Eigen::Matrix2d xhat = x.get_matrix();
        Eigen::Matrix2d Xhat = Eigen::Matrix2d::Zero(2,2);

        const double a = xhat(0,0);
        const double b = xhat(0,1);
        const double c = xhat(1,0);
        const double d = xhat(1,1);

        const double delta2 = std::pow(a - d, 2.0) + 4*b*c;
        const double m = std::exp((a + d)/2.0);

        if (delta2 < -1.0*eps)
        {
            // Ref [1] Corollary 2.4 Case 3.
            const double delta = 0.5*std::sqrt(std::abs(delta2));
            Xhat(0,0) = m*(std::cos(delta) + (a - d)/2.0*std::sin(delta)/delta);
            Xhat(0,1) = m*(b*std::sin(delta)/delta);
            Xhat(1,0) = m*(c*std::sin(delta)/delta);
            Xhat(1,1) = m*(std::cos(delta) - (a - d)/2.0*std::sin(delta)/delta);
        }
        else if (delta2 > 1.0*eps)
        {
            // Ref [1] Corollary 2.4 Case 2.
            const double delta = 0.5*std::sqrt(delta2);
            Xhat(0,0) = m*(std::cosh(delta) + (a - d)/2.0*std::sinh(delta)/delta);
            Xhat(0,1) = m*(b*std::sinh(delta)/delta);
            Xhat(1,0) = m*(c*std::sinh(delta)/delta);
            Xhat(1,1) = m*(std::cosh(delta) - (a - d)/2.0*std::sinh(delta)/delta);
        }
        else
        {
            // Ref [1] Corollary 2.4 Case 1.
            Xhat(0,0) = m*(1.0 + (a - d)/2.0);
            Xhat(0,1) = m*(b);
            Xhat(1,0) = m*(c);
            Xhat(1,1) = m*(1.0 - (a - d)/2.0);
        }

        return Lielab::domain::GLR(Xhat);
    }

    return exp_numerical(x);
}

template <>
Lielab::domain::RN exp(const Lielab::domain::rn& x)
{
    /*! \f{equation*}{ (\mathfrak{rn}) \rightarrow RN \f}
    
    Exponential function overload for \f$\mathfrak{rn}\f$.

    \f$\mathfrak{rn}\f$ is Abelian so the data is copied directly.
    
    Arguments
    ---------
    @param[in] la An instance of rn
    @param[out] out An instance of RN

    */

    return Lielab::domain::RN::from_vector(x.point);
}

template <>
Lielab::domain::SE exp(const Lielab::domain::se& x)
{
    /*! \f{equation*}{ (\mathfrak{se}) \rightarrow SE \f}
    
    Exponential function overload for \f$\mathfrak{se}\f$.

    For \f$y \in \mathfrak{se}(2)\f$ (shape 3), this uses [1]:

    \f{equation*}{y = \begin{bmatrix}
    w & x \\
    \mathbf{0} & 0
    \end{bmatrix}, \; \theta = \Vert w \Vert \f}

    \f{equation*}{\exp(y) = \begin{bmatrix}
    \exp(w) & Vx \\
    \mathbf{0} & 1
    \end{bmatrix} \f}

    where

    \f{equation*}{V = \frac{1}{\theta}\begin{bmatrix}
    \sin(\theta) & \cos(\theta) - 1 \\
    1 - \cos(\theta) & \sin(\theta)
    \end{bmatrix} \f}
    
    For \f$y \in \mathfrak{se}(3)\f$ (shape 4), this uses [1,2]:

    \f{equation*}{y = \begin{bmatrix}
    w & x \\
    \mathbf{0} & 0
    \end{bmatrix}, \; \theta = \Vert w \Vert \f}

    \f{equation*}{\exp(y) = \begin{bmatrix}
    \exp(w) & Vx \\
    \mathbf{0} & 1
    \end{bmatrix} \f}

    where

    \f{equation*}{V = \mathbf{I} + \frac{1 - \cos(\theta)}{\theta^2}\hat{w} + \frac{\theta - \sin(\theta)}{\theta^3}\hat{w}^2 \f}
    
    For \f$y \in \mathfrak{se}(4+)\f$ (shape 5 and greater), this uses the numerical procedure.

    Arguments
    ---------
    @param[in] x An instance of se
    @param[out] out An instance of SE

    References
    ----------

    [1] Ethan Eade. Lie groups for 2d and 3d transformations. Technical report,
        May 2017. [Online]. Available: http://ethaneade.com/lie.pdf.
    
    [2] Jose Luis Blanco-Claraco. A tutorial on se(3) transformation
        parameterizations and on-manifold optimization. arXiv preprint
        arXiv:2103.15980, 2021.

    */

    using Lielab::domain::RN;
    using Lielab::domain::SO;
    using Lielab::domain::SE;
    using Lielab::domain::rn;
    using Lielab::domain::so;
    using Lielab::domain::se;

    const size_t shape = x.get_shape();

    if (shape == 3)
    {
        // Source: Eade
        Lielab::domain::SE out = Lielab::domain::SE::identity(3);
        const auto [rn_c, so_c] = x.get_point();
        const Eigen::VectorXd rn_cbar = rn_c.get_vector();
        const Eigen::VectorXd so_cbar = so_c.get_vector();

        // SO component. Re-use SO calculation.
        std::get<SO>(out.point) = Lielab::functions::exp<Lielab::domain::so>(so_c);
        
        // R component
        const double theta = so_cbar(0);
        double cA = std::sin(theta)/theta;
        double cB = (1.0 - std::cos(theta))/theta;
        if (std::abs(theta) <= 1.0e-14)
        {
            cA = 1.0;
            cB = 0.0;
        }

        std::get<RN>(out.point) = Lielab::domain::RN::from_vector({cA*rn_cbar(0) - cB*rn_cbar(1), cB*rn_cbar(0) + cA*rn_cbar(1)});

        return out;
    }

    if (shape == 4)
    {
        // Sources: Eade and Blanco-Claraco
        Lielab::domain::SE out = Lielab::domain::SE::identity(4);
        const auto [rn_c, so_c] = x.get_point();
        const Eigen::VectorXd rn_cbar = rn_c.get_vector();
        const Eigen::VectorXd so_cbar = so_c.get_vector();
        const Eigen::MatrixXd what = so_c.get_matrix();
        const Eigen::MatrixXd what2 = what*what;

        // SO component. Re-use SO calculation.
        const Lielab::domain::so w(what);
        std::get<SO>(out.point) = Lielab::functions::exp<Lielab::domain::so>(w);
        
        // R component
        const double theta = std::sqrt(std::pow(so_cbar(0), 2.0) + std::pow(so_cbar(1), 2.0) + std::pow(so_cbar(2), 2.0));
        const double theta2 = std::pow(theta, 2.0);
        const double theta3 = std::pow(theta, 3.0);
        const double stheta = std::sin(theta);
        const double ctheta = std::cos(theta);
        const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(3, 3);

        double c1 = (1.0 - ctheta)/theta2;
        double c2 = (theta - stheta)/theta3;
        if (std::abs(theta) <= 1e-14)
        {
            c1 = 0.0;
            c2 = 0.0;
        }
        const Eigen::MatrixXd V = Id + c1*what + c2*what2;
        const Eigen::VectorXd Xbar = V*rn_cbar;

        std::get<RN>(out.point).unserialize(Xbar);

        return out;
    }

    return exp_numerical(x);
}

template <>
Lielab::domain::SO exp(const Lielab::domain::so& x)
{
    /*! \f{equation*}{ (\mathfrak{so}) \rightarrow SO \f}
    
    Exponential function overload for \f$\mathfrak{so}\f$.

    For \f$x \in \mathfrak{so}(2)\f$, this uses:

    \f{equation*}{\theta = \Vert x \Vert\f}

    \f{equation*}{\exp(x) = \begin{bmatrix}
    \cos \theta & - \sin \theta \\
    \sin \theta & \cos \theta
    \end{bmatrix} \f}

    For \f$x \in \mathfrak{so}(3)\f$, this uses the Euler-Rodriguez formula [1]:

    \f{equation*}{\theta = \Vert x \Vert, \; v = \frac{\theta}{2} \f}

    \f{equation*}{\exp(x) = \mathbf{I} + \frac{\sin(\theta)}{\theta}\hat{x} + \frac{1}{2}\frac{\sin^2(v)}{v^2}\hat{x}^2 \f}
    
    For \f$x \in \mathfrak{so}(4+)\f$, this uses the numerical procedure.

    Arguments
    ---------
    @param[in] x An instance of so
    @param[out] out An instance of SO

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000

    */

    using Lielab::domain::SO;

    const size_t shape = x.get_shape();

    if (shape == 2)
    {
        const Eigen::VectorXd v = x.get_vector();
        const double st = std::sin(v(0));
        const double ct = std::cos(v(0));
        Eigen::MatrixXd outdata = Eigen::MatrixXd(2, 2);
        outdata(0, 0) = ct;
        outdata(0, 1) = -st;
        outdata(1, 0) = st;
        outdata(1, 1) = ct;
        return SO(outdata);
    }

    if (shape == 3)
    {
        // Euler-Rodrigues formula for so(3). Source: Iserles
        const Eigen::Vector3d xbar = x.get_vector();
        const Eigen::MatrixXd xhat = x.get_matrix();
        const Eigen::MatrixXd xhat2 = xhat*xhat;
        const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);

        const double theta = std::sqrt(std::pow(xbar(0), 2.0) + std::pow(xbar(1), 2.0) + std::pow(xbar(2), 2.0));
        const double v = theta/2.0;
        const double v2 = std::pow(v, 2.0);

        double c1 = std::sin(theta)/theta;
        double c2 = 1.0/2.0*std::pow(std::sin(v), 2.0)/(v2);

        if (std::abs(theta) <= 1e-14)
        {
            c1 = 1.0;
            c2 = 0.5;
        }

        return SO(I + c1*xhat + c2*xhat2);
    }

    return exp_numerical(x);
}

// sp

// su

// log

template <>
Lielab::domain::cn log(const Lielab::domain::CN& g)
{
    /*! \f{equation*}{ (CN) \rightarrow \mathfrak{cn} \f}
    
    Logarithm function overload for \f$CN\f$.

    \f$CN\f$ is Abelian so the data is copied directly.
    
    Arguments
    ---------
    @param[in] g An instance of CN
    @param[out] out An instance of cn

    */

    return Lielab::domain::cn::from_complex_vector(g.point);
}

template <>
Lielab::domain::rn log(const Lielab::domain::RN& g)
{
    /*! \f{equation*}{ (RN) \rightarrow \mathfrak{rn} \f}
    
    Logarithm function overload for \f$RN\f$.

    \f$RN\f$ is Abelian so the data is copied directly.
    
    Arguments
    ---------
    @param[in] G An instance of RN
    @param[out] out An instance of rn

    */

    return Lielab::domain::rn::from_vector(g.serialize());
}

template <>
Lielab::domain::so log(const Lielab::domain::SO& g)
{
    /*! \f{equation*}{ (SO) \rightarrow \mathfrak{so} \f}

    Logarithm function overload for \f$SO\f$.

    For \f$W \in SO(2)\f$, this uses:

    \f{equation*}{\theta = \text{atan2}(\hat{g}_{10}, \hat{g}_{00})\f}

    \f{equation*}{\log(g) = \begin{bmatrix}
    0 & - \theta \\
    \theta & 0
    \end{bmatrix} \f}

    For \f$g \in SO(3)\f$, this uses [3]:
    tbd...
    // \f{equation*}{\hat{y} = \frac{1}{2}(\hat{g} - \hat{g}^T), \; \theta = \Vert \hat{y} \Vert \f}

    // \f{equation*}{\text{log}(g) = \frac{\sin^{-1}(\theta)}{\theta}\hat{y} \f}
    
    For \f$W \in SO(4+)\f$, this uses the numerical procedure.
    
    Arguments
    ---------
    @param[in] g An instance of SO
    @param[out] out An instance of so

    References
    ----------
    [1] Kenth Engø. On the bch-formula in so(3). BIT Numerical Mathematics,
        41:629–632, 2001
    
    [2] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000
    
    [3] Bullo, Francesco, and Andrew D. Lewis. Geometric control of mechanical
        systems: modeling, analysis, and design for simple mechanical control
        systems. Vol. 49. Springer, 2019.
    
    */

    const size_t shape = g.get_shape();

    if (shape == 2)
    {
        // Sources: Refs 1 and 2.
        Lielab::domain::so out(2);

        const Eigen::MatrixXd ghat = g.get_matrix();
        Eigen::VectorXd v = Eigen::VectorXd::Zero(1);
        v(0) = std::atan2(ghat(1,0), ghat(0,0));

        out.set_vector(v);
        return out;
    }

    if (shape == 3)
    {
        // Source: Ref 3. Proposition 5.7.
        Lielab::domain::so out(3);
        
        const Eigen::MatrixXd ghat = g.get_matrix();
        const double theta = std::acos((ghat.trace() - 1)/2.0);
        
        double c = theta/(2.0*std::sin(theta));
        if (std::abs(theta) <= 1e-14)
        {
            c = 0.5;
        }

        out.point = c*(ghat - ghat.transpose());

        return out;
    }

    return log_numerical(g);
}

template <>
Lielab::domain::se log(const Lielab::domain::SE& g)
{
    /*! \f{equation*}{ (SE) \rightarrow \mathfrak{se} \f}
    
    Logarithm function overload for \f$SE\f$.

    For \f$Y \in SE(2)\f$ (shape 3), this uses [1]:

    \f{equation*}{g = \begin{bmatrix}
    W & X \\
    \mathbf{0} & 1
    \end{bmatrix}, \; \theta = \Vert \log(W) \Vert \f}

    \f{equation*}{\log(g) = \begin{bmatrix}
    \log(W) & V^{-1} X \\
    \mathbf{0} & 0
    \end{bmatrix} \f}

    where

    \f{equation*}{V = \frac{1}{\theta}\begin{bmatrix}
    \sin(\theta) & \cos(\theta) - 1 \\
    1 - \cos(\theta) & \sin(\theta)
    \end{bmatrix} \f}
    
    For \f$Y \in SE(3)\f$ (shape 4), this uses [1,2]:

    \f{equation*}{g = \begin{bmatrix}
    W & X \\
    \mathbf{0} & 1
    \end{bmatrix}, \; \theta = \Vert \log(W) \Vert \f}

    \f{equation*}{\exp(y) = \begin{bmatrix}
    \log(W) & V^{-1}X \\
    \mathbf{0} & 0
    \end{bmatrix} \f}

    where

    \f{equation*}{V = \mathbf{I} + \frac{1 - \cos(\theta)}{\theta^2}\hat{w} + \frac{\theta - \sin(\theta)}{\theta^3}\hat{w}^2 \f}
    
    For \f$Y \in SE(4+)\f$ (shape 5 and greater), this uses the numerical procedure.

    Arguments
    ---------
    @param[in] g An instance of SE
    @param[out] out An instance of se
    
    References
    ----------

    [1] Ethan Eade. Lie groups for 2d and 3d transformations. Technical report,
        May 2017. [Online]. Available: http://ethaneade.com/lie.pdf.
    
    [2] Jose Luis Blanco-Claraco. A tutorial on se(3) transformation
        parameterizations and on-manifold optimization. arXiv preprint
        arXiv:2103.15980, 2021.
    */

    using Lielab::domain::RN;
    using Lielab::domain::SO;
    using Lielab::domain::SE;
    using Lielab::domain::rn;
    using Lielab::domain::so;
    using Lielab::domain::se;

    const size_t shape = g.get_shape();

    if (shape == 3)
    {
        // Source: Eade
        se out = se::zero(3);
        const auto [RN_c, SO_c] = g.get_point();

        // SO component. Re-use SO calculation.
        std::get<so>(out.point) = log<SO>(SO_c);
        
        // R component
        const Eigen::VectorXd RN_cbar = RN_c.serialize();
        const double theta = std::get<so>(out.point).get_vector()(0);
        const double x = RN_cbar(0);
        const double y = RN_cbar(1);

        double cA = std::sin(theta)/theta;
        double cB = (1.0 - std::cos(theta))/theta;
        
        if (std::abs(theta) <= 1e-14)
        {
            cA = 1.0;
            cB = 0.0;
        }
        
        const double den = std::pow(cA, 2.0) + std::pow(cB, 2.0);
        const double cAden = cA/den;
        const double cBden = cB/den;
        std::get<rn>(out.point).set_vector({cAden*x + cBden*y, -cBden*x + cAden*y});

        return out;
    }

    if (shape == 4)
    {
        // Sources: Eade and Blanco-Claraco
        se out = se::zero(4);
        const auto [RN_c, SO_c] = g.get_point();

        // so component. Re-use so calculation.
        std::get<so>(out.point) = log<SO>(SO_c);
        
        // R component
        const Eigen::VectorXd xbar = RN_c.serialize();
        const Eigen::VectorXd wbar = std::get<so>(out.point).get_vector();
        const Eigen::MatrixXd what = std::get<so>(out.point).get_matrix();
        const double theta = std::sqrt(std::pow(wbar(0), 2.0) + std::pow(wbar(1), 2.0) + std::pow(wbar(2), 2.0));
        const double theta2 = std::pow(theta, 2.0);
        const double theta3 = std::pow(theta, 3.0);
        const double stheta = std::sin(theta);
        const double ctheta = std::cos(theta);

        double cA = (1.0 - ctheta)/theta2;
        double cB = (theta - stheta)/theta3;

        if (theta <= 1e-14)
        {
            cA = 0.5;
            cB = 1.0/6.0;
        }

        const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(3, 3);
        const Eigen::MatrixXd V = Id + cA*what + cB*what*what;
        std::get<rn>(out.point).set_vector(V.inverse()*xbar);

        return out;
    }

    return log_numerical(g);
}

template <>
Lielab::domain::glr dexp(const Lielab::domain::glr& x, const int order)
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
Lielab::domain::glr dexp(const Lielab::domain::glr& x, const Lielab::domain::glr& y, const int order)
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
    
    const int shape = x.get_shape();

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
Lielab::domain::glr dexp(const Lielab::domain::se& y, const int order)
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
            // aw = 1.0;
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
Lielab::domain::se dexp(const Lielab::domain::se& x, const Lielab::domain::se& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
    
    Overloaded derivative of the exponential function for se.

    For \f$x \in \mathfrak{se}(2)\f$ (shape 3), this uses [1]:

    \f{equation*}{\text{dexp}_{x}(y) = \text{dexp}(x)\bar{y} \f}

    For \f$x \in \mathfrak{se}(3)\f$ (shape 4), this uses [1]:

    \f{equation*}{\text{dexp}_{x}(y) = \text{dexp}(x)\bar{y} \f}

    For \f$x \in \mathfrak{se}(4+)\f$ (shape 5+), this uses the numerical procedure.

    Arguments
    ---------
    @param[in] x First instance of se.
    @param[in] y Second instance of se.
    @param[in] order Order of the series expansion. Default 5.
    @param[out] out An instance of se.

    References
    ----------
    [1] Ethan Eade. Derivative of the exponential map. Technical report,
        November 2018. [Online]. Available: https://ethaneade.com/exp_diff.pdf.

    */

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dexp: Shapes of a and b must be equal.");
    }

    if (shape == 3)
    {
        const Lielab::domain::glr left = dexp(x);
        const Eigen::MatrixXd ybar = y.get_vector();
        return Lielab::domain::se::from_vector(left.get_matrix()*ybar);
    }

    if (shape == 4)
    {
        const Lielab::domain::glr left = dexp(x);
        const Eigen::MatrixXd ybar = y.get_vector();
        return Lielab::domain::se::from_vector(left.get_matrix()*ybar);
    }

    return dexp_numerical<Lielab::domain::se>(x, y, order);
}

template <>
Lielab::domain::glr dexp(const Lielab::domain::so& x, const int order)
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
    [2] Derived it myself.

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
        if (std::abs(theta) <= 1e-14)
        {
            c1 = 0.5;
        }

        // Source: Derived this myself. 3-point Lagrange polynomial with Lobatto node spacing fit at the points 0.0, 0.005, and 0.01.
        // Only interpolate starting from the midpoint, 0.005, for continuity but includes 0.01 for more accurate trend information.
        const double eps = std::numeric_limits<double>::epsilon();
        double c2 = (theta - std::sin(theta))/(theta3);
        if (theta < 10.0*eps)
        {
            c2 = 1.0/6.0;
        }
        else if (theta < 0.005 - 10.0*eps)
        {
            const double num = 20000.0*(1.0/6.0)/theta - 40000.0*0.16666645833579571/(theta - 0.005) + 20000.0*0.166665833335744/(theta - 0.01);
            const double den = 20000.0/theta - 40000.0/(theta - 0.005) + 20000.0/(theta - 0.01);
            c2 = num/den;
        }

        return Lielab::domain::glr(I + c1*ad1 + c2*ad2);
    }

    return dexp_numerical<Lielab::domain::so>(x, order);
}

template <>
Lielab::domain::so dexp(const Lielab::domain::so& x, const Lielab::domain::so& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{so} \f}
    
    Overloaded derivative of the exponential function for so.

    For \f$x \in \mathfrak{so}(2)\f$, this uses [1]:

    \f{equation*}{\text{dexp}_{x}(y) = y \f}

    For \f$x \in \mathfrak{so}(3)\f$, this uses [1]:

    \f{equation*}{\text{dexp}_{x}(y) = \text{dexp}(x)\bar{y} \f}

    For \f$a \in \mathfrak{so}(4+)\f$, this uses the numerical procedure.

    Arguments
    ---------
    @param[in] x First instance of so
    @param[in] y Second instance of so
    @param[in] order Order of the series expansion.
    @param[out] out An instance of so

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.

    */

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dexp: Shapes of x and y must be equal.");
    }

    if (shape == 2)
    {
        return y;
    }

    if (shape == 3)
    {
        const Lielab::domain::glr left = dexp(x);
        const Eigen::MatrixXd ybar = y.get_vector();
        Lielab::domain::so out(3);
        out.set_vector(left.get_matrix()*ybar);
        return out;
    }

    return dexp_numerical<Lielab::domain::so>(x, y, order);
}

template <>
Lielab::domain::glr dexp(const Lielab::domain::su& x, const int order)
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
Lielab::domain::su dexp(const Lielab::domain::su& x, const Lielab::domain::su& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    
    Overloaded derivative of the exponential function for su.

    For \f$x \in \mathfrak{su}(2)\f$, this uses:

    \f{equation*}{\text{dexp}_{x}(y) = \text{dexp}(x)\bar{y} \f}

    For \f$a \in \mathfrak{su}(3+)\f$ and \f$\mathfrak{su}(2)\f$, this uses the numerical procedure.

    Arguments
    ---------
    @param[in] x First instance of so
    @param[in] y Second instance of so
    @param[in] order Order of the series expansion.
    @param[out] out An instance of so

    References
    ----------
    Derived it myself - Mike Sparapany

    */

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dexp: Shapes of x and y must be equal.");
    }

    // if (shape == 1)
    // {
    //     return b;
    // }

    if (shape == 2)
    {
        const Lielab::domain::glr left = dexp(x);
        const Eigen::MatrixXd ybar = y.get_vector();
        Lielab::domain::su out(2);
        out.set_vector(left.get_matrix()*ybar);
        return out;
    }

    return dexp_numerical<Lielab::domain::su>(x, y, order);
}



template <>
Lielab::domain::glr dexpinv(const Lielab::domain::glr& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overloaded inverse derivative of the exponential function for glr.

    Arguments
    ---------
    @param[in] x First instance of glr
    @param[in] order Order of the series expansion.
    @param[out] out An instance of glr

    References
    ----------
    [1] Derived it myself.

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

        const double c1 = -0.5;
        double c2;
        if (std::abs(theta) <= 1e-14)
        {
            c2 = -1.0/12.0;
        }
        else if (theta <= -1e-14)
        {
            c2 = (std::sqrt(-theta)*std::cos(std::sqrt(-theta)/2.0)/std::sin(std::sqrt(-theta)/2.0) - 2.0)/(2.0*theta);
        }
        else
        {
            c2 = (std::sqrt(theta)*std::cosh(std::sqrt(theta)/2.0)/std::sinh(std::sqrt(theta)/2.0) - 2.0)/(2.0*theta);
        }

        return Lielab::domain::glr(I + c1*ad1 + c2*ad2);
    }

    return dexpinv_numerical<Lielab::domain::glr>(x, order);
}

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::glr& x, const Lielab::domain::glr& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overloaded inverse derivative of the exponential function for glr.

    Arguments
    ---------
    @param[in] x First instance of glr
    @param[in] y Second instance of glr
    @param[in] order Order of the series expansion.
    @param[out] out An instance of glr

    References
    ----------
    [1] Derived it myself.

    */

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dexpinv: Shapes of x and y must be equal.");
    }

    if (shape == 2)
    {
        const Eigen::MatrixXd adxhat = dexpinv<Lielab::domain::glr>(x).get_matrix();
        const Eigen::MatrixXd yhat = y.get_vector();
        Lielab::domain::glr out(2);
        out.set_vector(adxhat*yhat);
        return out;
    }

    return dexpinv_numerical<Lielab::domain::glr>(x, y, order);
}

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::se& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
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
        const Lielab::domain::glr left_inv = Lielab::functions::dexp(y);
        const Eigen::MatrixXd A = left_inv.get_matrix().block(0, 0, 2, 2);
        const Eigen::VectorXd v = left_inv.get_matrix().block(0, 2, 2, 1);
        const Eigen::MatrixXd Ainv = A.inverse();

        Eigen::MatrixXd left = Eigen::MatrixXd::Zero(3, 3);
        left.block(0, 0, 2, 2) = Ainv;
        left.block(0, 2, 2, 1) = -Ainv*v;
        left(2, 2) = 1.0;
        return Lielab::domain::glr(left);
    }

    if (shape == 4)
    {
        // Source: Eade
        const Lielab::domain::glr left_inv = Lielab::functions::dexp(y);
        const Eigen::MatrixXd A = left_inv.get_matrix().block(0, 0, 3, 3);
        const Eigen::MatrixXd B = left_inv.get_matrix().block(0, 3, 3, 3);
        const Eigen::MatrixXd Ainv = A.inverse();

        Eigen::MatrixXd left = Eigen::MatrixXd::Zero(6, 6);

        left.block(0, 0, 3, 3) = Ainv;
        left.block(0, 3, 3, 3) = -Ainv*B*Ainv;
        left.block(3, 3, 3, 3) = Ainv;
        return Lielab::domain::glr(left);
    }

    return dexpinv_numerical<Lielab::domain::se>(y, order);
}

template <>
Lielab::domain::se dexpinv(const Lielab::domain::se& x, const Lielab::domain::se& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
    
    Overloaded derivative of the exponential function for se.

    For \f$a \in \mathfrak{se}(2)\f$ (shape 3), this uses [1]:

    \f{equation*}{\text{dexp}^{-1}_{x}(y) = \text{dexp}^{-1}(x)\bar{y} \f}

    For \f$a \in \mathfrak{se}(3)\f$ (shape 4), this uses [1]:

    \f{equation*}{\text{dexp}^{-1}_{x}(y) = \text{dexp}^{-1}(x)\bar{y} \f}

    For \f$a \in \mathfrak{se}(4+)\f$ (shape 5+), this uses the numerical procedure.

    Arguments
    ---------
    @param[in] x First instance of se
    @param[in] y Second instance of se
    @param[in] order Order of the series expansion.
    @param[out] out An instance of se

    References
    ----------
    [1] Ethan Eade. Derivative of the exponential map. Technical report,
        November 2018. [Online]. Available: https://ethaneade.com/exp_diff.pdf.

    */

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dexpinv: Shapes of x and y must be equal.");
    }

    if (shape == 3)
    {
        const Lielab::domain::glr left = dexpinv(x);
        const Eigen::MatrixXd ybar = y.get_vector();
        return Lielab::domain::se::from_vector(left.get_matrix()*ybar);
    }

    if (shape == 4)
    {
        const Lielab::domain::glr left = dexpinv(x);
        const Eigen::MatrixXd ybar = y.get_vector();
        return Lielab::domain::se::from_vector(left.get_matrix()*ybar);
    }

    return dexpinv_numerical<Lielab::domain::se>(x, y, order);
}

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::so& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overloaded inverse derivative of the exponential function for so.

    For \f$x \in \mathfrak{so}(3)\f$, this uses [1]:

    \f{equation*}{\theta = \Vert x \Vert, \; v = \frac{\theta}{2} \f}

    \f{equation*}{\text{dexp}^{-1}(x) = \mathbf{I} - \frac{1}{2}\hat{x} - \frac{\theta \cot(v) - 2}{2 \theta^2}\hat{x}^2 \f}

    Arguments
    ---------
    @param[in] x First instance of so
    @param[in] order Order of the series expansion.
    @param[out] out An instance of glr

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.
    [2] Derived some limits myself.

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
        const Eigen::MatrixXd xhat = x.get_matrix();
        const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);

        const double theta = std::sqrt(std::pow(xbar(0), 2.0) + std::pow(xbar(1), 2.0) + std::pow(xbar(2), 2.0));
        const double v = theta/2.0;
        const double sv = std::sin(v);
        const double cv = std::cos(v);
        const double theta2 = std::pow(theta, 2.0);
        const Eigen::MatrixXd xhat2 = xhat*xhat;

        // Source: Derived this myself. 3-point Lagrange polynomial with Lobatto node spacing fit at the points 0.0, 0.005, and 0.01.
        // Only interpolate starting from the midpoint, 0.005, for continuity but includes 0.01 for more accurate trend information.
        const double eps = std::numeric_limits<double>::epsilon();
        double c2 = (theta*cv/sv - 2.0)/(2.0*theta2);
        if (theta < 10.0*eps)
        {
            c2 = -1.0/12.0;
        }
        else if (theta < 0.005 - 10.0*eps)
        {
            const double num = 20000.0*(-1.0/12.0)/theta - 40000.0*-0.08333336805499414/(theta - 0.005) + 20000.0*-0.08333347222166942/(theta - 0.01);
            const double den = 20000.0/theta - 40000.0/(theta - 0.005) + 20000.0/(theta - 0.01);
            c2 = num/den;
        }

        const Eigen::MatrixXd left = I - 1.0/2.0*xhat - c2*xhat2;
        return Lielab::domain::glr(left);
    }

    return dexpinv_numerical<Lielab::domain::so>(x, order);
}

template <>
Lielab::domain::so dexpinv(const Lielab::domain::so& x, const Lielab::domain::so& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{so} \f}
    
    Overloaded inverse derivative of the exponential function for so.

    For \f$x \in \mathfrak{so}(3)\f$, this uses [1]:

    \f{equation*}{\text{dexp}^{-1}_{x}(y) = \text{dexp}^{-1}(x)\bar{y} \f}

    For \f$a \in \mathfrak{so}(4+)\f$ and \f$\mathfrak{so}(2)\f$, this uses the numerical procedure.

    Arguments
    ---------
    @param[in] x First instance of so
    @param[in] y Second instance of so
    @param[in] order Order of the series expansion.
    @param[out] out An instance of so

    References
    ----------
    [1] Arieh Iserles, Hans Z Munthe-Kaas, Syvert P Nørsett, and Antonella
        Zanna. Lie-group methods. Acta numerica, 9:215–365, 2000.

    */

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dexpinv: Shapes of x and y must be equal.");
    }

    if (shape == 2)
    {
        return y;
    }

    if (shape == 3)
    {
        const Lielab::domain::glr left = dexpinv(x);
        const Eigen::MatrixXd ybar = y.get_vector();
        Lielab::domain::so out(3);
        out.set_vector(left.get_matrix()*ybar);
        return out;
    }

    return dexpinv_numerical<Lielab::domain::so>(x, y, order);
}

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::su& x, const int order)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overloaded inverse derivative of the exponential function for so.

    For \f$x \in \mathfrak{su}(2)\f$, this uses:

    // \f{equation*}{\theta = \Vert x \Vert, \; v = \frac{\theta}{2} \f}

    // \f{equation*}{\text{dexp}^{-1}(x) = \mathbf{I} - \frac{1}{2}\hat{x} - \frac{\theta \cot(v) - 2}{2 \theta^2}\hat{x}^2 \f}

    Arguments
    ---------
    @param[in] x First instance of su
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

        const double theta = std::sqrt(std::pow(xbar(0), 2.0) + std::pow(xbar(1), 2.0) + std::pow(xbar(2), 2.0));
        const double st = std::sin(theta);
        const double ct = std::cos(theta);
        const double theta2 = std::pow(theta, 2.0);

        double c2 = (2.0*theta*ct/st - 2.0)/(8.0*theta2);
        if (std::abs(theta) <= 1e-14)
        {
            c2 = -1.0/12.0;
        }

        const Eigen::MatrixXd left = I - 1.0/2.0*ad1 - c2*ad2;
        return Lielab::domain::glr(left);
    }

    return dexpinv_numerical<Lielab::domain::su>(x, order);
}

template <>
Lielab::domain::su dexpinv(const Lielab::domain::su& x, const Lielab::domain::su& y, const int order)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    
    Overloaded inverse derivative of the exponential function for su.

    For \f$x \in \mathfrak{su}(2)\f$, this uses:

    // \f{equation*}{\text{dexp}^{-1}_{x}(y) = \text{dexp}^{-1}(x)\bar{y} \f}

    For \f$a \in \mathfrak{su}(3+)\f$ this uses the numerical procedure.

    Arguments
    ---------
    @param[in] x First instance of su
    @param[in] y Second instance of su
    @param[in] order Order of the series expansion.
    @param[out] out An instance of su

    References
    ----------
    Derived it myself - Mike Sparapany

    */

    const int shape = x.get_shape();

    if (shape != y.get_shape())
    {
        throw Lielab::utils::InputError("dexpinv: Shapes of x and y must be equal.");
    }

    // if (shape == 1)
    // {
    //     return y;
    // }

    if (shape == 2)
    {
        const Lielab::domain::glr left = dexpinv(x);
        const Eigen::MatrixXd ybar = y.get_vector();
        Lielab::domain::su out(2);
        out.set_vector(left.get_matrix()*ybar);
        return out;
    }

    return dexpinv_numerical<Lielab::domain::su>(x, y, order);
}

// Composite overloads

template <>
Lielab::domain::CompositeGroup exp_numerical(const Lielab::domain::CompositeAlgebra& x)
{
    /*!
    * CompositeAlgebra exponential overload.
    */
    
    using namespace Lielab::domain;

    CompositeGroup out;

    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(exp_numerical(_element));
        }, element);
    }

    return out;
}

template <>
Lielab::domain::CompositeGroup exp(const Lielab::domain::CompositeAlgebra& x)
{
    /*!
    * CompositeAlgebra exponential overload.
    */
    
    using namespace Lielab::domain;

    CompositeGroup out;

    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(exp(_element));
        }, element);
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra log_numerical(const Lielab::domain::CompositeGroup& g)
{
    /*!
     * CompositeGroup logarithm overload
     */

    using namespace Lielab::domain;

    CompositeAlgebra out;
    
    for (const auto& element : g.point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(log_numerical(_element));
        }, element);
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra log(const Lielab::domain::CompositeGroup& g)
{
    /*!
     * CompositeGroup logarithm overload
     */

    using namespace Lielab::domain;

    CompositeAlgebra out;
    
    for (const auto& element : g.point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(log(_element));
        }, element);
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dexp_numerical(const Lielab::domain::CompositeAlgebra& x, const int order)
{
    /*!
    * CompositeAlgebra derivative of the exponential numerical overload.
    */

    using namespace Lielab::domain;

    CompositeAlgebra out;

    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(dexp_numerical(_element, order));
        }, element);
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dexp(const Lielab::domain::CompositeAlgebra& x, const int order)
{
    /*!
    * CompositeAlgebra derivative of the exponential overload.
    */

    using namespace Lielab::domain;

    CompositeAlgebra out;

    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(dexp(_element, order));
        }, element);
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dexp_numerical(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order)
{
    /*!
    * CompositeAlgebra derivative of the exponential numerical overload.
    */

    using Lielab::domain::CompositeAlgebra;
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(x, y), "Unable to take numerical dexp of topologically inconsistent algebras: (" + x.to_string() + ") !≅ (" + y.to_string() + ").");

    CompositeAlgebra out;
    int index = 0;
    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            out.point.push_back(dexp_numerical(_element, std::get<other_t>(y.point[index]), order));
        }, element);
        index++;
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dexp(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order)
{
    /*!
    * CompositeAlgebra derivative of the exponential overload.
    */

    using Lielab::domain::CompositeAlgebra;
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(x, y), "Unable to take dexp of topologically inconsistent algebras: (" + x.to_string() + ") !≅ (" + y.to_string() + ").");

    CompositeAlgebra out;
    int index = 0;
    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            out.point.push_back(dexp(_element, std::get<other_t>(y.point[index]), order));
        }, element);
        index++;
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dexpinv_numerical(const Lielab::domain::CompositeAlgebra& x, const int order)
{
    /*!
    * CompositeAlgebra dexpinv_numerical overload
    */

    using namespace Lielab::domain;

    CompositeAlgebra out;

    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(dexpinv_numerical(_element, order));
        }, element);
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dexpinv(const Lielab::domain::CompositeAlgebra& x, const int order)
{
    /*!
    * CompositeAlgebra dexpinv overload
    */

    using namespace Lielab::domain;

    CompositeAlgebra out;

    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(dexpinv(_element, order));
        }, element);
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dexpinv_numerical(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order)
{
    /*!
    * CompositeAlgebra dexpinv_numerical overload
    */

    using Lielab::domain::CompositeAlgebra;
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(x, y), "Unable to take numerical dexpinv of topologically inconsistent algebras: (" + x.to_string() + ") !≅ (" + y.to_string() + ").");

    CompositeAlgebra out;
    int index = 0;
    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            out.point.push_back(dexpinv_numerical(_element, std::get<other_t>(y.point[index]), order));
        }, element);
        index++;
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dexpinv(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order)
{
    /*!
    * CompositeAlgebra dexpinv overload
    */

    using Lielab::domain::CompositeAlgebra;
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(x, y), "Unable to take dexpinv of topologically inconsistent algebras: (" + x.to_string() + ") !≅ (" + y.to_string() + ").");

    CompositeAlgebra out;
    int index = 0;
    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            out.point.push_back(dexpinv(_element, std::get<other_t>(y.point[index]), order));
        }, element);
        index++;
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dlog_numerical(const Lielab::domain::CompositeAlgebra& x, const int order)
{
    return dexpinv_numerical<Lielab::domain::CompositeAlgebra, Lielab::domain::CompositeAlgebra>(x, order);
}

template <>
Lielab::domain::CompositeAlgebra dlog(const Lielab::domain::CompositeAlgebra& x, const int order)
{
    return dexpinv<Lielab::domain::CompositeAlgebra, Lielab::domain::CompositeAlgebra>(x, order);
}

template <>
Lielab::domain::CompositeAlgebra dlog_numerical(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order)
{
    return dexpinv_numerical<Lielab::domain::CompositeAlgebra>(x, y, order);
}

template <>
Lielab::domain::CompositeAlgebra dlog(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order)
{
    return dexpinv<Lielab::domain::CompositeAlgebra>(x, y, order);
}

template <>
Lielab::domain::CompositeAlgebra dloginv_numerical(const Lielab::domain::CompositeAlgebra& x, const int order)
{
    return dexp_numerical<Lielab::domain::CompositeAlgebra, Lielab::domain::CompositeAlgebra>(x, order);
}

template <>
Lielab::domain::CompositeAlgebra dloginv(const Lielab::domain::CompositeAlgebra& x, const int order)
{
    return dexp<Lielab::domain::CompositeAlgebra, Lielab::domain::CompositeAlgebra>(x, order);
}

template <>
Lielab::domain::CompositeAlgebra dloginv_numerical(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order)
{
    return dexp_numerical<Lielab::domain::CompositeAlgebra>(x, y, order);
}

template <>
Lielab::domain::CompositeAlgebra dloginv(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y, const int order)
{
    return dexp<Lielab::domain::CompositeAlgebra>(x, y, order);
}

}
