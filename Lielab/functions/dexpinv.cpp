#include "dexpinv.hpp"

#include "dexp.hpp"
#include "littlead.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::glr & x, const size_t order)
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
Lielab::domain::glr dexpinv(const Lielab::domain::glr & x, const Lielab::domain::glr & y, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overloaded inverse derivative of the exponential function for glr.

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
Lielab::domain::glr dexpinv(const Lielab::domain::so & x, const size_t order)
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

        double c2 = (theta*cv/sv - 2.0)/(2.0*theta2);
        if (std::abs(theta) <= 1e-14)
        {
            c2 = -1.0/12.0;
        }

        const Eigen::MatrixXd left = I - 1.0/2.0*xhat - c2*xhat2;
        return Lielab::domain::glr(left);
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

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("dexpinv: Shapes of a and b must be equal.");
    }

    if (shape == 2)
    {
        return b;
    }

    if (shape == 3)
    {
        const Lielab::domain::glr left = dexpinv(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::so out(3);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    return dexpinv_numerical<Lielab::domain::so>(a, b, order);
}

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::se & y, const size_t order)
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

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("dexpinv: Shapes of a and b must be equal.");
    }

    if (shape == 3)
    {
        const Lielab::domain::glr left = dexpinv(a);
        const Eigen::MatrixXd y = b.get_vector();
        return Lielab::domain::se::from_vector(left.get_matrix()*y);
    }

    if (shape == 4)
    {
        const Lielab::domain::glr left = dexpinv(a);
        const Eigen::MatrixXd y = b.get_vector();
        return Lielab::domain::se::from_vector(left.get_matrix()*y);
    }

    return dexp_numerical<Lielab::domain::se>(a, b, order);
}

template <>
Lielab::domain::glr dexpinv(const Lielab::domain::su & x, const size_t order)
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
Lielab::domain::su dexpinv(const Lielab::domain::su & a, const Lielab::domain::su & b, const size_t order)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{su} \f}
    
    Overloaded inverse derivative of the exponential function for su.

    For \f$x \in \mathfrak{su}(2)\f$, this uses:

    // \f{equation*}{\text{dexp}^{-1}_{x}(y) = \text{dexp}^{-1}(x)\bar{y} \f}

    For \f$a \in \mathfrak{su}(3+)\f$ this uses the numerical procedure.

    Arguments
    ---------
    @param[in] a First instance of su
    @param[in] b Second instance of su
    @param[in] order Order of the series expansion.
    @param[out] out An instance of su

    References
    ----------
    Derived it myself - Mike Sparapany

    */

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("dexpinv: Shapes of a and b must be equal.");
    }

    // if (shape == 1)
    // {
    //     return b;
    // }

    if (shape == 2)
    {
        const Lielab::domain::glr left = dexpinv(a);
        const Eigen::MatrixXd y = b.get_vector();
        Lielab::domain::su out(2);
        out.set_vector(left.get_matrix()*y);
        return out;
    }

    return dexpinv_numerical<Lielab::domain::su>(a, b, order);
}

template <>
Lielab::domain::CompositeAlgebra dexpinv(const Lielab::domain::CompositeAlgebra & a, const Lielab::domain::CompositeAlgebra & b, const size_t order)
{
    /*!
    * CompositeAlgebra dexpinv overload
    */

    Lielab::domain::CompositeAlgebra out;

    for (size_t ii = 0; ii < a.space.size(); ii++)
    {
        const size_t ind = a.space[ii].index();
        if (ind == Lielab::domain::CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::cn>(a.space[ii]),
                                                           std::get<Lielab::domain::cn>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::glr>(a.space[ii]),
                                                           std::get<Lielab::domain::glr>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::glc>(a.space[ii]),
                                                           std::get<Lielab::domain::glc>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::rn>(a.space[ii]),
                                                           std::get<Lielab::domain::rn>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::se>(a.space[ii]),
                                                           std::get<Lielab::domain::se>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::so>(a.space[ii]),
                                                           std::get<Lielab::domain::so>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::sp>(a.space[ii]),
                                                           std::get<Lielab::domain::sp>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::su>(a.space[ii]),
                                                           std::get<Lielab::domain::su>(b.space[ii]), order));
        }
    }

    return out;
}

}
