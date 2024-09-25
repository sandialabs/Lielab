#ifndef _LIELAB_FUNCTIONS_EXP_HPP
#define _LIELAB_FUNCTIONS_EXP_HPP

#include "../domain.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
namespace functions
{

template<typename LA>
Lielab::domain::lieiii<LA> exp_numerical(const LA & la)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
    
    This is the main exponential function. Computes the exponential in a
    computationally intensive numerical process.
    
    Arguments
    ---------
    @param[in] la An instance of g
    @param[out] out An instance of G

    */

    return (la.get_matrix()).exp();
}

template<typename LA>
Lielab::domain::lieiii<LA> exp(const LA & la)
{
    /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
    
    Catch-all function for the exponential function. Will always use the
    numerical procedure.
    
    Arguments
    ---------
    @param[in] la An instance of g
    @param[out] out An instance of G

    */
    
    return exp_numerical(la);
}

template<>
Lielab::domain::CN exp(const Lielab::domain::cn & la)
{
    /*! \f{equation*}{ (\mathfrak{cn}) \rightarrow CN \f}
    
    Exponential function overload for \f$\mathfrak{cn}\f$.

    \f$\mathfrak{cn}\f$ is Abelian so the data is copied directly.
    
    Arguments
    ---------
    @param[in] la An instance of cn
    @param[out] out An instance of CN

    */

    return la._data;
}

template<>
Lielab::domain::RN exp(const Lielab::domain::rn & la)
{
    /*! \f{equation*}{ (\mathfrak{rn}) \rightarrow RN \f}
    
    Exponential function overload for \f$\mathfrak{rn}\f$.

    \f$\mathfrak{rn}\f$ is Abelian so the data is copied directly.
    
    Arguments
    ---------
    @param[in] la An instance of rn
    @param[out] out An instance of RN

    */

    return la._data;
}

template<>
Lielab::domain::SO exp(const Lielab::domain::so & x)
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

    static const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(3,3);
    static const Eigen::MatrixXd ZZ = Eigen::MatrixXd::Zero(3,3);

    if (x.shape == 2)
    {
        const Eigen::VectorXd v = x.get_vector();
        const double st = std::sin(v(0));
        const double ct = std::cos(v(0));
        Eigen::MatrixXd outdata = Eigen::MatrixXd::Zero(2,2);
        outdata(0,0) = ct;
        outdata(0,1) = -st;
        outdata(1,0) = st;
        outdata(1,1) = ct;
        return Lielab::domain::SO(outdata);
    }

    if (x.shape == 3)
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

        const Eigen::MatrixXd data = I + c1*xhat + c2*xhat2;
        return Lielab::domain::SO(data);
    }

    return exp_numerical(x);
}

template<>
Lielab::domain::SE exp(const Lielab::domain::se & y)
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
    @param[in] y An instance of se
    @param[out] out An instance of SE

    References
    ----------

    [1] Ethan Eade. Lie groups for 2d and 3d transformations. Technical report,
        May 2017. [Online]. Available: http://ethaneade.com/lie.pdf.
    
    [2] Jose Luis Blanco-Claraco. A tutorial on se(3) transformation
        parameterizations and on-manifold optimization. arXiv preprint
        arXiv:2103.15980, 2021.

    */

    if (y.shape == 3)
    {
        // Source: Eade
        Lielab::domain::SE out(3);
        const Eigen::VectorXd ybar = y.get_vector();
        const Eigen::MatrixXd yhat = y.get_matrix();
        const Eigen::MatrixXd what = yhat.block(0, 0, 2, 2);

        // SO component. Re-use SO calculation.
        const Lielab::domain::so w(what);
        const Lielab::domain::SO W = Lielab::functions::exp<Lielab::domain::so>(w);
        out._data.block(0, 0, 2, 2) = W.get_matrix();
        
        // R component
        const double theta = ybar(2);
        double cA = std::sin(theta)/theta;
        double cB = (1.0 - std::cos(theta))/theta;
        if (std::abs(theta) <= 1.0e-14)
        {
            cA = 1.0;
            cB = 0.0;
        }
        out._data(0, 2) = cA*ybar(0) - cB*ybar(1);
        out._data(1, 2) = cB*ybar(0) + cA*ybar(1);

        return out;
    }

    if (y.shape == 4)
    {
        // Sources: Eade and Blanco-Claraco
        Lielab::domain::SE out(4);
        const Eigen::VectorXd ybar = y.get_vector();
        const Eigen::MatrixXd yhat = y.get_matrix();
        const Eigen::MatrixXd what = yhat.block(0, 0, 3, 3);
        const Eigen::MatrixXd what2 = what*what;

        // SO component. Re-use SO calculation.
        const Lielab::domain::so w(what);
        const Lielab::domain::SO W = Lielab::functions::exp<Lielab::domain::so>(w);
        out._data.block(0, 0, 3, 3) = W.get_matrix();
        
        // R component
        const Eigen::VectorXd xbar = ybar(Eigen::seqN(0, 3));
        const double theta = std::sqrt(std::pow(ybar(3), 2.0) + std::pow(ybar(4), 2.0) + std::pow(ybar(5), 2.0));
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
        const Eigen::VectorXd Xbar = V*xbar;
        out._data(0, 3) = Xbar(0);
        out._data(1, 3) = Xbar(1);
        out._data(2, 3) = Xbar(2);

        return out;
    }

    return exp_numerical(y);
}


}
}

#endif
