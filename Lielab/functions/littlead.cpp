#include "littlead.hpp"

#include "commutator.hpp"

#include "Lielab/utils.hpp"

#include <vector>

namespace Lielab::functions
{

template <>
Lielab::domain::glr ad(const Lielab::domain::cn & a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Adjoint function overload for cn.
    
    This returns:

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^j_a = \mathbf{0}, \; j \neq 0 \f}

    Arguments
    ---------
    @param[in] a First element.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of g.
    
    */

    const size_t dim = a.get_dimension();

    if (p == 0)
    {
        return Lielab::domain::glr(Eigen::MatrixXd::Identity(dim, dim));
    }

    return Lielab::domain::glr(dim);
}

template <>
Lielab::domain::cn ad(const Lielab::domain::cn & a, const Lielab::domain::cn & b, const int p)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{cn} \f}
    
    Adjoint function overload for cn.
    
    This returns:

    \f{equation*}{\text{ad}^0_a b = b \f}

    \f{equation*}{\text{ad}^j_a b = \mathbf{0}, \; j \neq 0 \f}

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
        throw Lielab::utils::InputError("ad: Shapes of a and b must be equal.");
    }

    if (p == 0)
    {
        return b;
    }

    return std::complex<double>(0.0, 0.0)*b;
}

template <>
Lielab::domain::glr ad(const Lielab::domain::glr & x, const int p)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overload for the little adjoint function of glr.

    Arguments
    ---------
    @param[in] a Element of glr.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of glr.

    References
    ----------
    [1] I made it up.
    
    */
    
    const size_t shape = x.get_shape();

    if (shape == 2)
    {
        if (p == 0)
        {
            const Eigen::MatrixXd adpxhat = Eigen::MatrixXd::Identity(4,4);
            return Lielab::domain::glr(adpxhat);
        }

        const Eigen::Matrix2d xhat = x.get_matrix();
        const double a = xhat(0,0);
        const double b = xhat(0,1);
        const double c = xhat(1,0);
        const double d = xhat(1,1);
        const double theta = std::pow(a - d, 2.0) + 4.0*b*c;

        Eigen::MatrixXd ad1 = Eigen::MatrixXd::Zero(4, 4);

        ad1(0,1) = -c;
        ad1(0,2) = b;
        ad1(1,0) = -b;
        ad1(1,1) = a - d;
        ad1(1,3) = b;
        ad1(2,0) = c;
        ad1(2,2) = d - a;
        ad1(2,3) = -c;
        ad1(3,1) = c;
        ad1(3,2) = -b;

        const double n = static_cast<double>((p-1)/2);
        const double coeff = std::pow(theta, n);

        if ((p%2) == 1)
        {
            return Lielab::domain::glr(coeff*ad1);
        }
        else
        {
            return Lielab::domain::glr(coeff*ad1*ad1);
        }
    }

    return ad_numerical<Lielab::domain::glr>(x, p);
}

template <>
Lielab::domain::glr ad(const Lielab::domain::glr & a, const Lielab::domain::glr & b, const int p)
{
    /*! \f{equation*}{ (\mathfrak{glr}, \mathfrak{glr}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overload for the little adjoint function of glr.

    For \f$a \in \mathfrak{glr}(2) \f$, this uses a known formula.

    Otherwise, this uses a numerical formula.

    Arguments
    ---------
    @param[in] a First element of glr.
    @param[in] b Second element of glr.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of glr.
    
    */

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("ad: Shapes of a and b must be equal.");
    }

    if (shape == 2)
    {
        const Eigen::MatrixXd adjaphat = ad<Lielab::domain::glr>(a, p).get_matrix();
        const Eigen::VectorXd bbar = b.get_vector();
        Lielab::domain::glr out(shape);
        out.set_vector(adjaphat*bbar);
        return out;
    }
    
    return ad_numerical<Lielab::domain::glr>(a, b, p);
}

// template <>
// Lielab::domain::glr ad(const Lielab::domain::glc & x, const int p)
// {
//     /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
//     Overload for the little adjoint function of glc.

//     Arguments
//     ---------
//     @param[in] a Element of glc.
//     @param[in] p Power of the adjoint.
//     @param[out] out An instance of glr.

//     References
//     ----------
//     [1] I made it up.
    
//     */
    
//     if (x.shape == 2)
//     {
//         if (p == 0)
//         {
//             const Eigen::MatrixXd adpxhat = Eigen::MatrixXd::Identity(8, 8);
//             return Lielab::domain::glr(adpxhat);
//         }

//         const Eigen::Matrix2cd xhat = x.get_matrix();
//         const std::complex<double> a = xhat(0,0);
//         const std::complex<double> b = xhat(0,1);
//         const std::complex<double> c = xhat(1,0);
//         const std::complex<double> d = xhat(1,1);
//         const std::complex<double> theta = std::pow(a - d, 2.0) + 4.0*b*c;

//         Eigen::MatrixXd ad1 = Eigen::MatrixXd::Zero(8, 8);

//         ad1(0,2) = -c.real();
//         ad1(0,3) = c.imag();
//         ad1(0,4) = b.real();
//         ad1(0,5) = -b.imag();
//         ad1(1,2) = -c.imag();
//         ad1(1,3) = -c.real();
//         ad1(1,4) = b.imag();
//         ad1(1,5) = b.real();
//         ad1(2,0) = -b.real();
//         ad1(2,1) = b.imag();
//         ad1(2,2) = a.real() - d.real();
//         ad1(2,3) = -a.imag() + d.imag();
//         ad1(2,6) = b.real();
//         ad1(2,7) = -b.imag();
//         ad1(3,0) = -b.imag();
//         ad1(3,1) = -b.real();
//         ad1(3,2) = a.imag() - d.imag();
//         ad1(3,3) = a.real() - d.real();
//         ad1(3,6) = b.imag();
//         ad1(3,7) = b.real();
//         ad1(4,0) = c.real();
//         ad1(4,1) = -c.imag();
//         ad1(4,4) = -a.real() + d.real();
//         ad1(4,5) = a.imag() - d.imag();
//         ad1(4,6) = -c.real();
//         ad1(4,7) = c.imag();
//         ad1(5,0) = c.imag();
//         ad1(5,1) = c.real();
//         ad1(5,4) = -a.imag() + d.imag();
//         ad1(5,5) = -a.real() + d.real();
//         ad1(5,6) = -c.imag();
//         ad1(5,7) = -c.real();
//         ad1(6,2) = c.real();
//         ad1(6,3) = -c.imag();
//         ad1(6,4) = -b.real();
//         ad1(6,5) = b.imag();
//         ad1(7,2) = c.imag();
//         ad1(7,3) = c.real();
//         ad1(7,4) = -b.imag();
//         ad1(7,5) = -b.real();

//         // ad1 = ([[ Z, A, B, Z],
//         //         [-B, C, Z, B],
//         //         [-A, Z,-C, A],
//         //         [ Z,-A,-B, Z]])

//         // ad3 = ([[0, -4*A**2*B + A*C**2, -4*A*B**2 + B*C**2, 0],
//         //         [2*A*B**2 - B*(-2*A*B + C**2), -2*A*B*C + C*(-2*A*B + C**2), 0, -2*A*B**2 + B*(-2*A*B + C**2)],
//         //         [2*A**2*B - A*(-2*A*B + C**2), 0, 2*A*B*C - C*(-2*A*B + C**2), -2*A**2*B + A*(-2*A*B + C**2)],
//         //         [0, 4*A**2*B - A*C**2, 4*A*B**2 - B*C**2, 0]])

//         // TODO: the rest of this function
//         // const double n = static_cast<double>((p-1)/2);
//         // const double coeff = 1.0;  //std::pow(theta, n); // TODO:

//         // if ((p%2) == 1)
//         // {
//         //     return Lielab::domain::glr(coeff*ad1);
//         // }
//         // else
//         // {
//         //     return Lielab::domain::glr(coeff*ad1*ad1);
//         // }
//     }

//     return ad_numerical<Lielab::domain::glc>(x, p);
// }

// template <>
// Lielab::domain::glc ad(const Lielab::domain::glc & a, const Lielab::domain::glc & b, const int p)
// {
//     /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}, \mathbb{R}) \rightarrow \mathfrak{glc} \f}
    
//     Overload for the little adjoint function of glr.

//     For \f$a \in \mathfrak{glc}(2) \f$, this uses a known formula.

//     Otherwise, this uses a numerical formula.

//     Arguments
//     ---------
//     @param[in] a First element of glc.
//     @param[in] b Second element of glc.
//     @param[in] p Power of the adjoint.
//     @param[out] out An instance of glc.
    
//     */

//     if (a.shape != b.shape)
//     {
//         throw Lielab::utils::InputError("ad: Shapes of a and b must be equal.");
//     }

//     const ptrdiff_t shape = a.shape;
//     if (shape == 2)
//     {
//         const Eigen::MatrixXd adjaphat = ad<Lielab::domain::glc>(a, p).get_matrix();
//         const Eigen::VectorXd bbar = b.get_vector();
//         Lielab::domain::glc out(shape);
//         out.set_vector(adjaphat*bbar);
//         return out;
//     }
    
//     return ad_numerical<Lielab::domain::glc>(a, b, p);
// }

template <>
Lielab::domain::glr ad(const Lielab::domain::rn & a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Adjoint function overload for rn.
    
    This returns:

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^j_a = \mathbf{0}, \; j \neq 0 \f}

    Arguments
    ---------
    @param[in] a First element.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of g.
    
    */

    const size_t dim = a.get_dimension();

    if (p == 0)
    {
        return Lielab::domain::glr(Eigen::MatrixXd::Identity(dim, dim));
    }

    return Lielab::domain::glr(dim);
}

template <>
Lielab::domain::rn ad(const Lielab::domain::rn & a, const Lielab::domain::rn & b, const int p)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
    
    Adjoint function overload for rn.
    
    This returns:

    \f{equation*}{\text{ad}^0_a b = b \f}

    \f{equation*}{\text{ad}^j_a b = \mathbf{0}, \; j \neq 0 \f}

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
        throw Lielab::utils::InputError("ad: Shapes of a and b must be equal.");
    }

    if (p == 0)
    {
        return b;
    }

    return 0.0*b;
}

template <>
Lielab::domain::glr ad(const Lielab::domain::so & a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overload for the little adjoint function of so.

    For \f$a \in \mathfrak{so}(2)\f$, this uses:

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^j_a = \mathbf{0}, j \neq 0 \f}

    For \f$a \in \mathfrak{so}(3)\f$, this uses [1]:

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^1_a = \hat{a} \f}

    \f{equation*}{\text{ad}^2_a = \hat{a}^2 \f}

    \f{equation*}{\text{ad}^3_a = -\Vert a \Vert^2 \text{ad}^1_a \f}

    Arguments
    ---------
    @param[in] a Element of so.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of glr.

    References
    ----------
    [1] Ethan Eade. Derivative of the exponential map. Technical report,
        November 2018. [Online]. Available: https://ethaneade.com/exp_diff.pdf.

    TODO
    ----
        - Change recursive call for p > 3 to use mod+remainder division.
          This will prevent freezing up when the stack grows large.
    
    */

    const size_t shape = a.get_shape();

    if (shape == 2)
    {
        if (p == 0)
        {
            return Lielab::domain::glr(Eigen::MatrixXd::Identity(1, 1));
        }

        return Lielab::domain::glr(1);
    }

    if (shape == 3)
    {
        Lielab::domain::glr out(3);

        if (p == 0)
        {
            out.data(0, 0) = 1.0;
            out.data(1, 1) = 1.0;
            out.data(2, 2) = 1.0;
            return out;
        }

        const Eigen::VectorXd abar = a.get_vector();
        const double wmag2 = std::pow(abar(0), 2.0) + std::pow(abar(1), 2.0) + std::pow(abar(2), 2.0);
        const double n = static_cast<double>((p-1)/2);
        const double coeff = std::pow(-wmag2, n);

        const Eigen::MatrixXd ad1 = a.get_matrix();

        if ((p%2) == 1)
        {
            out.data = coeff*ad1;
        }
        else
        {
            out.data = coeff*ad1*ad1;
        }

        return out;
    }

    return ad_numerical<Lielab::domain::so>(a, p);
}

template <>
Lielab::domain::so ad(const Lielab::domain::so & a, const Lielab::domain::so & b, const int p)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overload for the little adjoint function of so.

    For \f$a \in \mathfrak{so}(2) \vee \mathfrak{so}(3) \f$, this uses a known formula.

    Otherwise, this uses a numerical formula.

    Arguments
    ---------
    @param[in] a First element of so.
    @param[in] b Second element of so.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of so.
    
    */

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("ad: Shapes of a and b must be equal.");
    }

    if (p == 0)
    {
        return b;
    }

    if (shape == 2 || shape == 3)
    {
        const Eigen::MatrixXd adjaphat = ad<Lielab::domain::so>(a, p).get_matrix();
        const Eigen::VectorXd bbar = b.get_vector();
        Lielab::domain::so out(shape);
        out.set_vector(adjaphat*bbar);
        return out;
    }
    
    return ad_numerical<Lielab::domain::so>(a, b, p);
}

template <>
Lielab::domain::glr ad(const Lielab::domain::su & a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overload for the little adjoint function of su.

    For \f$a \in \mathfrak{su}(2)\f$, this uses:

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^1_a = \begin{bmatrix}
    0 & -2 a_3 & a_2 \\
    2 a_3 & 0 & -2 a_1 \\
    -2 a_2 & 2 x_1 & 0
    \end{bmatrix} \f}

    \f{equation*}{\text{ad}^2_a = (\text{ad}_a)^2 \f}

    \f{equation*}{\text{ad}^3_a = -4 \Vert a \Vert^2 \text{ad}^1_a \f}

    Arguments
    ---------
    @param[in] a Element of su.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of glr.

    References
    ----------
    Derived this myself - Mike Sparapany
    
    */

    const size_t shape = a.get_shape();

    if (shape == 2)
    {
        Lielab::domain::glr out(3);

        if (p == 0)
        {
            out.data(0, 0) = 1.0;
            out.data(1, 1) = 1.0;
            out.data(2, 2) = 1.0;
            return out;
        }

        const Eigen::VectorXd abar = a.get_vector();
        const double wmag2 = std::pow(abar(0), 2.0) + std::pow(abar(1), 2.0) + std::pow(abar(2), 2.0);
        const double n = static_cast<double>((p-1)/2);
        const double coeff = std::pow(-4.0*wmag2, n);

        Eigen::MatrixXd ad1 = Eigen::MatrixXd::Zero(3, 3);
        ad1(0,1) = -2.0*abar(2);
        ad1(0,2) =  2.0*abar(1);
        ad1(1,2) = -2.0*abar(0);
        ad1(1,0) =  2.0*abar(2);
        ad1(2,0) = -2.0*abar(1);
        ad1(2,1) =  2.0*abar(0);

        if ((p%2) == 1)
        {
            out.data = coeff*ad1;
        }
        else
        {
            out.data = coeff*ad1*ad1;
        }

        return out;
    }

    return ad_numerical<Lielab::domain::su>(a, p);
}

template <>
Lielab::domain::su ad(const Lielab::domain::su & a, const Lielab::domain::su & b, const int p)
{
    /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overload for the little adjoint function of su.

    For \f$a \in \mathfrak{su}(2) \f$, this uses a known formula.

    Otherwise, this uses a numerical formula.

    Arguments
    ---------
    @param[in] a First element of su.
    @param[in] b Second element of su.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of su.
    
    */

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("ad: Shapes of a and b must be equal.");
    }

    if (p == 0)
    {
        return b;
    }

    if (shape == 2)
    {
        const Eigen::MatrixXd adjaphat = ad<Lielab::domain::su>(a, p).get_matrix();
        const Eigen::VectorXd bbar = b.get_vector();
        Lielab::domain::su out(shape);
        out.set_vector(adjaphat*bbar);
        return out;
    }
    
    return ad_numerical<Lielab::domain::su>(a, b, p);
}

template <>
Lielab::domain::glr ad(const Lielab::domain::se & a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overload for the little adjoint function of se.

    For \f$a \in \mathfrak{se}(2)\f$ (shape 3), this uses (Eq. 80 from Eade is incorrect (use theta^2 instead of theta^3)) [1]:

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^1_a = \begin{bmatrix}
    0 & -\theta & y \\
    \theta & 0 & -x \\
    0 & 0 & 0
    \end{bmatrix} \f}

    \f{equation*}{\text{ad}^2_a = \begin{bmatrix}
    -\theta^2 & 0 & x \theta \\
    0 & -\theta^2 & y \theta \\
    0 & 0 & 0
    \end{bmatrix} \f}

    \f{equation*}{\text{ad}^3_a = - \theta^2 \text{ad}^1_a \f}

    For \f$a \in \mathfrak{se}(3)\f$ (shape 4), this uses [1]:

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^1_a = \begin{bmatrix}
    \hat{w} & \hat{x} \\
    \mathbf{0} & \hat{w}
    \end{bmatrix} \f}

    \f{equation*}{\text{ad}^2_a = \begin{bmatrix}
    \hat{w}^2 & \hat{w}\hat{x} + \hat{x}\hat{w} \\
    \mathbf{0} & \hat{w}^2
    \end{bmatrix} \f}

    \f{equation*}{\text{ad}^4_a = -\theta^2 \text{ad}^2_a - 2 (\bar{w}^T \bar{x}) \begin{bmatrix}
    \mathbf{0} & \hat{w} \\
    \mathbf{0} & \mathbf{0}
    \end{bmatrix} \text{ad}^1_a \f}

    Arguments
    ---------
    @param[in] a Element of se.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of glr.

    References
    ----------
    [1] Ethan Eade. Derivative of the exponential map. Technical report,
        November 2018. [Online]. Available: https://ethaneade.com/exp_diff.pdf.
    
    Also derived myself - Mike Sparapany

    TODO
    ----
        - Change recursive call for p > 3 to use mod+remainder division.
          This will prevent freezing up when the stack grows large.

    */

    const size_t shape = a.get_shape();

    if (shape == 3)
    {
        Lielab::domain::glr out(3);

        if (p == 0)
        {
            out.data(0, 0) = 1.0;
            out.data(1, 1) = 1.0;
            out.data(2, 2) = 1.0;
            return out;
        }

        const Eigen::VectorXd abar = a.get_vector();
        const double x = abar(0);
        const double y = abar(1);
        const double theta = abar(2);
        const double theta2 = std::pow(theta, 2.0);
        const double n = static_cast<double>((p-1)/2);
        const double coeff = std::pow(-theta2, n);

        if ((p%2) == 1)
        {
            out.data(0, 1) = -theta;
            out.data(1, 0) = theta;
            out.data(0, 2) = y;
            out.data(1, 2) = -x;
        }
        else
        {
            out.data(0, 0) = -theta2;
            out.data(1, 1) = -theta2;
            out.data(0, 2) = theta*x;
            out.data(1, 2) = theta*y;
        }
        
        out.data = coeff*out.data;
        return out;
    }

    if (shape == 4)
    {
        Lielab::domain::glr out(6);

        if (p == 0)
        {
            out.data(0, 0) = 1.0;
            out.data(1, 1) = 1.0;
            out.data(2, 2) = 1.0;
            out.data(3, 3) = 1.0;
            out.data(4, 4) = 1.0;
            out.data(5, 5) = 1.0;
            return out;
        }

        const Eigen::VectorXd abar = a.get_vector();
        const Eigen::MatrixXd ahat = a.get_matrix();
        const Eigen::MatrixXd what = ahat.block(0, 0, 3, 3);
        Lielab::domain::so x(3);
        x.set_vector(ahat.block(0, 3, 3, 1));
        const Eigen::MatrixXd xhat = x.get_matrix();
        Eigen::MatrixXd ad1 = Eigen::MatrixXd::Zero(6,6);
        ad1.block(0, 0, 3, 3) = what;
        ad1.block(0, 3, 3, 3) = xhat;
        ad1.block(3, 3, 3, 3) = what;

        if (p == 1)
        {
            out.data = ad1;
            return out;
        }

        const Eigen::MatrixXd what2 = what*what;
        Eigen::MatrixXd ad2 = Eigen::MatrixXd::Zero(6,6);
        ad2.block(0, 0, 3, 3) = what2;
        ad2.block(0, 3, 3, 3) = what*xhat + xhat*what;
        ad2.block(3, 3, 3, 3) = what2;

        if (p == 2)
        {
            out.data = ad2;
            return out;
        }

        // p >= 3
        const double theta2 = std::pow(abar(3), 2.0) + std::pow(abar(4), 2.0) + std::pow(abar(5), 2.0);
        const double innerprod = abar(0)*abar(3) + abar(1)*abar(4) + abar(2)*abar(5);
        Eigen::MatrixXd upperright = Eigen::MatrixXd::Zero(6, 6);
        upperright.block(0, 3, 3, 3) = what;

        Eigen::MatrixXd adm3 = Eigen::MatrixXd::Identity(6, 6);
        Eigen::MatrixXd adm2 = ad1;
        Eigen::MatrixXd adm1 = ad2;
        Eigen::MatrixXd adcurrent = -theta2*adm2 - 2.0*innerprod*upperright*adm3;

        ptrdiff_t ii = 3;
        while (ii < p)
        {
            adm3 = adm2;
            adm2 = adm1;
            adm1 = adcurrent;
            adcurrent = -theta2*adm2 - 2.0*innerprod*upperright*adm3;
            ii += 1;
        }

        out.data = adcurrent;
        return out;
    }

    return ad_numerical<Lielab::domain::se>(a, p);
}

template <>
Lielab::domain::se ad(const Lielab::domain::se & a, const Lielab::domain::se & b, const int p)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{glr} \f}
    
    Overload for the little adjoint function of se.

    For \f$a \in \mathfrak{se}(2) \vee \mathfrak{se}(3) \f$ (shapes 3 and 4), this uses a known formula.

    Otherwise, this uses a numerical formula.

    Arguments
    ---------
    @param[in] a First element of se.
    @param[in] b Second element of se.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of se.
    
    */

    const size_t shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("ad: Shapes of a and b must be equal.");
    }

    if (p == 0)
    {
        return b;
    }

    if (shape == 3 || shape == 4)
    {
        const Eigen::MatrixXd adjaphat = ad<Lielab::domain::se>(a, p).get_matrix();
        const Eigen::VectorXd bbar = b.get_vector();
        Lielab::domain::se out(shape);
        out.set_vector(adjaphat*bbar);
        return out;
    }
    
    return ad_numerical<Lielab::domain::se>(a, b, p);
}

}
