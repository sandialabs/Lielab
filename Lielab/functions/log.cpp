#include "log.hpp"

#include "Lielab/domain.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::functions
{

template <>
Lielab::domain::cn log(const Lielab::domain::CN& G)
{
    /*! \f{equation*}{ (CN) \rightarrow \mathfrak{cn} \f}
    
    Logarithm function overload for \f$CN\f$.

    \f$CN\f$ is Abelian so the data is copied directly.
    
    Arguments
    ---------
    @param[in] G An instance of CN
    @param[out] out An instance of cn

    */

    return Lielab::domain::cn::from_complex_vector(G.data);
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
Lielab::domain::so log(const Lielab::domain::SO& W)
{
    /*! \f{equation*}{ (SO) \rightarrow \mathfrak{so} \f}

    Logarithm function overload for \f$SO\f$.

    For \f$W \in SO(2)\f$, this uses:

    \f{equation*}{\theta = \text{atan2}(\hat{W}_{10}, \hat{W}_{00})\f}

    \f{equation*}{\log(W) = \begin{bmatrix}
    0 & - \theta \\
    \theta & 0
    \end{bmatrix} \f}

    For \f$W \in SO(3)\f$, this uses [1, 2]:

    \f{equation*}{\hat{y} = \frac{1}{2}(\hat{W} - \hat{W}^T), \; \theta = \Vert \hat{y} \Vert \f}

    \f{equation*}{\text{log}(W) = \frac{\sin^{-1}(\theta)}{\theta}\hat{y} \f}
    
    For \f$W \in SO(4+)\f$, this uses the numerical procedure.
    
    Arguments
    ---------
    @param[in] W An instance of SO
    @param[out] out An instance of so

    Notes
    -----
    Engø's formula is a local solution. Do not use when theta >= pi/2

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

    const size_t shape = W.get_shape();

    if (shape == 2)
    {
        // Sources: Refs 1 and 2.
        Lielab::domain::so out(2);

        const Eigen::MatrixXd What = W.get_matrix();
        Eigen::VectorXd v = Eigen::VectorXd::Zero(1);
        v(0) = std::atan2(What(1,0), What(0,0));

        out.set_vector(v);
        return out;
    }

    if (shape == 3)
    {
        // Source: Ref 3. Proposition 5.7.
        Lielab::domain::so out(3);
        
        const Eigen::MatrixXd What = W.get_matrix();
        const double theta = std::acos((What.trace() - 1)/2.0);
        
        double c = theta/(2.0*std::sin(theta));
        if (std::abs(theta) <= 1e-14)
        {
            c = 0.5;
        }

        out.data = c*(What - What.transpose());

        return out;
    }

    return log_numerical(W);
}

template<>
Lielab::domain::se log(const Lielab::domain::SE& Y)
{
    /*! \f{equation*}{ (SE) \rightarrow \mathfrak{se} \f}
    
    Logarithm function overload for \f$SE\f$.

    For \f$Y \in SE(2)\f$ (shape 3), this uses [1]:

    \f{equation*}{Y = \begin{bmatrix}
    W & X \\
    \mathbf{0} & 1
    \end{bmatrix}, \; \theta = \Vert \log(W) \Vert \f}

    \f{equation*}{\log(Y) = \begin{bmatrix}
    \log(W) & V^{-1} X \\
    \mathbf{0} & 0
    \end{bmatrix} \f}

    where

    \f{equation*}{V = \frac{1}{\theta}\begin{bmatrix}
    \sin(\theta) & \cos(\theta) - 1 \\
    1 - \cos(\theta) & \sin(\theta)
    \end{bmatrix} \f}
    
    For \f$Y \in SE(3)\f$ (shape 4), this uses [1,2]:

    \f{equation*}{Y = \begin{bmatrix}
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
    @param[in] Y An instance of SE
    @param[out] out An instance of se
    
    References
    ----------

    [1] Ethan Eade. Lie groups for 2d and 3d transformations. Technical report,
        May 2017. [Online]. Available: http://ethaneade.com/lie.pdf.
    
    [2] Jose Luis Blanco-Claraco. A tutorial on se(3) transformation
        parameterizations and on-manifold optimization. arXiv preprint
        arXiv:2103.15980, 2021.
    */

    const size_t shape = Y.get_shape();

    if (shape == 3)
    {
        // Source: Eade
        Lielab::domain::se out = Lielab::domain::se::from_shape(3);
        const Eigen::MatrixXd Yhat = Y.get_matrix();
        const Eigen::MatrixXd What = Yhat.block(0, 0, 2, 2);

        // SO component. Re-use SO calculation.
        const Lielab::domain::SO W(What);
        const Lielab::domain::so w = Lielab::functions::log<Lielab::domain::SO>(W);
        out.data.block(0, 0, 2, 2) = w.get_matrix();
        
        // R component
        const Eigen::VectorXd wbar = w.get_vector();
        const double x = Yhat(0, 2);
        const double y = Yhat(1, 2);
        const double theta = wbar(0);

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
        out.data(0, 2) =  cAden*x + cBden*y;
        out.data(1, 2) = -cBden*x + cAden*y;

        return out;
    }

    if (shape == 4)
    {
        // Sources: Eade and Blanco-Claraco
        Lielab::domain::se out = Lielab::domain::se::from_shape(4);
        const Eigen::MatrixXd Yhat = Y.get_matrix();
        const Eigen::MatrixXd What = Yhat.block(0, 0, 3, 3);

        // so component. Re-use so calculation.
        const Lielab::domain::SO W(What);
        const Lielab::domain::so w = Lielab::functions::log<Lielab::domain::SO>(W);
        out.data.block(0, 0, 3, 3) = w.get_matrix();
        
        // R component
        const Eigen::VectorXd wbar = w.get_vector();
        const Eigen::MatrixXd what = w.get_matrix();
        const Eigen::Vector3d Xbar(Yhat(0, 3), Yhat(1, 3), Yhat(2, 3));
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
        const Eigen::VectorXd xbar = V.inverse()*Xbar;
        out.data(0, 3) = xbar(0);
        out.data(1, 3) = xbar(1);
        out.data(2, 3) = xbar(2);

        return out;
    }

    return log_numerical(Y);
}

Lielab::domain::CompositeAlgebra log(const Lielab::domain::CompositeGroup& M)
{
    /*!
     * CompositeGroup logarithm overload
     */

    Lielab::domain::CompositeAlgebra out;
    
    for (size_t ii = 0; ii < M.space.size(); ii++)
    {
        const size_t ind = M.space[ii].index();
        if (ind == Lielab::domain::CompositeGroup::INDEX_CN)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::CN>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeGroup::INDEX_GLR)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::GLR>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeGroup::INDEX_GLC)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::GLC>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeGroup::INDEX_RN)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::RN>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeGroup::INDEX_SE)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::SE>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeGroup::INDEX_SO)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::SO>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeGroup::INDEX_SP)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::SP>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeGroup::INDEX_SU)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::SU>(M.space[ii])));
        }
    }

    return out;
}

}
