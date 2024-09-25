#ifndef LIELAB_FUNCTIONS_LITTLEAD_HPP_
#define LIELAB_FUNCTIONS_LITTLEAD_HPP_

#include <vector>

#include "../abstract.hpp"
#include "../utils.hpp"
#include "commutator.hpp"

namespace Lielab
{
namespace functions
{

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::gl ad_numerical(const LA & a, const int p = 1)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    Adjoint function on Lie algebras.
    
    This computes \f$\text{ad}_a^p\f$ in an entirely numerical manner. For basis
    \f$e_i \in \mathfrak{g}\f$, the Lie algebra structure constants are generated with:

    \f{equation*}{[e_i, e_j] = C_{ijk} e_k\f}

    Then the first adjoint of \f$a\f$ is found with:

    \f{equation*}{[\text{ad}_a]_{kj} = C_{ajk} \f}

    Powers of the adjoint are then simply

    \f{equation*}{\text{ad}^p_a = \Pi_p \text{ad}_a \f}

    Arguments
    ---------
    @param[in] a First element.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of g.

    Notes
    -----

    Has some known shortcuts for performance.

    1. Power 0 is identity.

    \f{equation*}{\text{ad}^0_a = \mathbf{I}, \, \forall a \in \mathfrak{g} \f}

    2. Abelian Lie algebras return 0, except for power 0.

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^j_a = \mathbf{0}, \, \text{if} \, \mathfrak{g} \in \mathfrak{Abelian}(\mathfrak{g}) \f}

    TODO
    ----
        - This function generates the structure constants with each call, thus
          making it unsuitable for heavy numerical use.
        - There exist formula for adjoint power multiples of 2. This could
          accelerate this procedure by quite a bit at higher powers.
    
    */

    const ptrdiff_t dim = a.get_dimension();
    
    // Shortcut for power 0 adjoints.
    if (p == 0)
    {
        return Lielab::domain::gl(Eigen::MatrixXd::Identity(dim, dim));
    }

    // Shortcut for Abelian Lie algebras.
    if (a.abelian)
    {
        return Lielab::domain::gl(dim);
    }

    const Eigen::VectorXd ahat = a.get_vector();
    const ptrdiff_t shape = a.shape;
    
    // Generate basis.
    std::vector<LA> basis;
    for (ptrdiff_t ii = 0; ii < dim; ii++)
    {
        basis.push_back(LA::basis(ii, shape));
    }

    // Generate structure constants.
    std::vector<std::vector<Eigen::VectorXd>> C;
    for (ptrdiff_t ii = 0; ii < dim; ii++)
    {
        std::vector<Eigen::VectorXd> vecj;
        for (ptrdiff_t jj = 0; jj < dim; jj++)
        {
            const Eigen::VectorXd veck = commutator(basis[ii], basis[jj]).get_vector();
            vecj.push_back(veck);
        }
        C.push_back(vecj);
    }

    // Generate first adjoint matrix.
    Eigen::MatrixXd adja = Eigen::MatrixXd::Zero(dim, dim);
    for (ptrdiff_t ii = 0; ii < dim; ii++)
    {
        for (ptrdiff_t jj = 0; jj < dim; jj++)
        {
            for (ptrdiff_t kk = 0; kk < dim; kk++)
            {
                adja(jj, ii) += ahat(kk)*C[kk][ii](jj);
            }
        }
    }

    // Raise adjoint matrix to the specified power.
    Eigen::MatrixXd adjap = Eigen::MatrixXd::Identity(dim, dim);
    for (ptrdiff_t ii = 0; ii < p; ii++)
    {
        adjap = adja*adjap;
    }

    return Lielab::domain::gl(adjap);
}

template <Lielab::abstract::LieAlgebra LA>
Lielab::domain::gl ad(const LA & a, const int p = 1)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    Catch-all function for the adjoint function on Lie algebras.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a First element.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of g.

    */
    
    return ad_numerical<LA>(a, p);
}

template <Lielab::abstract::LieAlgebra LA>
LA ad_numerical(const LA & a, const LA & b, const int p = 1)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Adjoint function on Lie algebras.
    
    This computes \f$\text{ad}_a^p b\f$ in an entirely numerical manner:

    \f{equation*}{\text{ad}^0_a b = b \f}

    \f{equation*}{\text{ad}^1_a b = [a, b] \f}

    \f{equation*}{\text{ad}^2_a b = [a, [a, b]] \f}

    ... and so on.

    Arguments
    ---------
    @param[in] a First element.
    @param[in] b Second element.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of g.
    
    */

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("ad_numerical: Shapes of a and b must be equal.");
    }

    LA adjb = b;

    for (ptrdiff_t ii = 0; ii < p; ii++)
    {
        adjb = commutator(a, adjb);
    }

    return adjb;
}

template <Lielab::abstract::LieAlgebra LA>
LA ad(const LA & a, const LA & b, const int p = 1)
{
    /*! \f{equation*}{ (\mathfrak{g}, \mathfrak{g}, \mathbb{R}) \rightarrow \mathfrak{g} \f}
    
    Catch-all function for the adjoint function on Lie algebras.
    Will always use the numerical procedure.

    Arguments
    ---------
    @param[in] a First element.
    @param[in] b Second element.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of g.

    */
    
    return ad_numerical<LA>(a, b, p);
}

template <>
Lielab::domain::gl ad(const Lielab::domain::cn & a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
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

    const ptrdiff_t dim = a.get_dimension();

    if (p == 0)
    {
        return Lielab::domain::gl(Eigen::MatrixXd::Identity(dim, dim));
    }

    return Lielab::domain::gl(dim);
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

    if (a.shape != b.shape)
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
Lielab::domain::gl ad(const Lielab::domain::rn & a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
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

    const ptrdiff_t dim = a.get_dimension();

    if (p == 0)
    {
        return Lielab::domain::gl(Eigen::MatrixXd::Identity(dim, dim));
    }

    return Lielab::domain::gl(dim);
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

    if (a.shape != b.shape)
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
Lielab::domain::gl ad(const Lielab::domain::so & a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    Overload for the little adjoint function of so.

    For \f$a \in \mathfrak{so}(2)\f$, this uses:

    \f{equation*}{\text{ad}^j_a = \mathbf{I} \f}

    For \f$a \in \mathfrak{so}(3)\f$, this uses [1]:

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^1_a = \hat{a} \f}

    \f{equation*}{\text{ad}^2_a = \hat{a}^2 \f}

    \f{equation*}{\text{ad}^3_a = -\Vert a \Vert^2 \text{ad}^1_a \f}

    Arguments
    ---------
    @param[in] a Element of so.
    @param[in] p Power of the adjoint.
    @param[out] out An instance of gl.

    References
    ----------
    [1] Ethan Eade. Derivative of the exponential map. Technical report,
        November 2018. [Online]. Available: https://ethaneade.com/exp_diff.pdf.

    TODO
    ----
        - Change recursive call for p > 3 to use mod+remainder division.
          This will prevent freezing up when the stack grows large.
    
    */

    if (a.shape == 2)
    {
        if (p == 0)
        {
            return Lielab::domain::gl(Eigen::MatrixXd::Identity(1, 1));
        }

        return Lielab::domain::gl(1);
    }

    if (a.shape == 3)
    {
        Lielab::domain::gl out(3);

        if (p == 0)
        {
            out(0,0) = 1.0;
            out(1,1) = 1.0;
            out(2,2) = 1.0;
            return out;
        }

        const Eigen::MatrixXd ahat = a.get_matrix();
        if (p == 1)
        {
            out._data = ahat;
            return out;
        }

        const Eigen::MatrixXd ahat2 = ahat*ahat;
        if (p == 2)
        {
            out._data = ahat2;
            return out;
        }

        const Eigen::VectorXd abar = a.get_vector();
        const double wmag2 = std::pow(abar(0), 2.0) + std::pow(abar(1), 2.0) + std::pow(abar(2), 2.0);
        return -wmag2*ad(a, p - 2);
    }

    return ad_numerical<Lielab::domain::so>(a, p);
}

template <>
Lielab::domain::so ad(const Lielab::domain::so & a, const Lielab::domain::so & b, const int p)
{
    /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
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

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("ad: Shapes of a and b must be equal.");
    }

    if (p == 0)
    {
        return b;
    }

    const ptrdiff_t shape = a.shape;
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
Lielab::domain::gl ad(const Lielab::domain::se & a, const int p)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
    Overload for the little adjoint function of se.

    For \f$a \in \mathfrak{se}(2)\f$ (shape 3), this uses [1]:

    \f{equation*}{\text{ad}^0_a = \mathbf{I} \f}

    \f{equation*}{\text{ad}^1_a = \begin{bmatrix}
    0 & -\theta & y \\
    \theta & 0 & -x \\
    0 & 0 & 0
    \end{bmatrix} \f}

    \f{equation*}{\text{ad}^2_a = \begin{bmatrix}
    0 & -\theta & y \\
    \theta & 0 & -x \\
    0 & 0 & 0
    \end{bmatrix} \f}

    \f{equation*}{\text{ad}^3_a = -\Vert \theta \Vert^3 \text{ad}^1_a \f}

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
    @param[out] out An instance of gl.

    References
    ----------
    [1] Ethan Eade. Derivative of the exponential map. Technical report,
        November 2018. [Online]. Available: https://ethaneade.com/exp_diff.pdf.

    TODO
    ----
        - Change recursive call for p > 3 to use mod+remainder division.
          This will prevent freezing up when the stack grows large.
        - se(3) (shape 4) is especially bad with higher powers since
          it uses 2 recursive calls each.

    */

    if (a.shape == 3)
    {
        Lielab::domain::gl out(3);

        if (p == 0)
        {
            out(0, 0) = 1.0;
            out(1, 1) = 1.0;
            out(2, 2) = 1.0;
            return out;
        }

        const Eigen::VectorXd abar = a.get_vector();
        const double x = abar(0);
        const double y = abar(1);
        const double theta = abar(2);
        if (p == 1)
        {
            out(0, 1) = -theta;
            out(1, 0) = theta;
            out(0, 2) = y;
            out(1, 2) = -x;
            return out;
        }

        if (p == 2)
        {
            const double theta2 = std::pow(theta, 2.0);
            out(0, 0) = -theta2;
            out(1, 1) = -theta2;
            out(0, 2) = theta*x;
            out(1, 2) = theta*y;
            return out;
        }
        
        const double theta3 = std::pow(theta, 3.0);
        return -theta3*ad(a, p - 2);
    }

    if (a.shape == 4)
    {
        Lielab::domain::gl out(6);

        if (p == 0)
        {
            out(0, 0) = 1.0;
            out(1, 1) = 1.0;
            out(2, 2) = 1.0;
            out(3, 3) = 1.0;
            out(4, 4) = 1.0;
            out(5, 5) = 1.0;
            return out;
        }

        const Eigen::VectorXd abar = a.get_vector();
        const Eigen::MatrixXd ahat = a.get_matrix();
        const Eigen::MatrixXd what = ahat.block(0, 0, 3, 3);
        Lielab::domain::so x(3);
        x.set_vector(ahat.block(0, 3, 3, 1));
        const Eigen::MatrixXd xhat = x.get_matrix();
        if (p == 1)
        {
            out._data.block(0, 0, 3, 3) = what;
            out._data.block(0, 3, 3, 3) = xhat;
            out._data.block(3, 3, 3, 3) = what;
            return out;
        }

        if (p == 2)
        {
            const Eigen::MatrixXd what2 = what*what;
            out._data.block(0, 0, 3, 3) = what2;
            out._data.block(0, 3, 3, 3) = what*xhat + xhat*what;
            out._data.block(3, 3, 3, 3) = what2;
            return out;
        }

        const double theta2 = std::pow(abar(3), 2.0) + std::pow(abar(4), 2.0) + std::pow(abar(5), 2.0);
        const double innerprod = abar(0)*abar(3) + abar(1)*abar(4) + abar(2)*abar(5);
        const Eigen::MatrixXd adlower2 = ad(a, p - 2).get_matrix();
        const Eigen::MatrixXd adlower3 = ad(a, p - 3).get_matrix();
        Eigen::MatrixXd upperright = Eigen::MatrixXd::Zero(6, 6);
        upperright.block(0, 3, 3, 3) = what;

        return -theta2*adlower2 - 2.0*innerprod*upperright*adlower3;
    }

    return ad_numerical<Lielab::domain::se>(a, p);
}

template <>
Lielab::domain::se ad(const Lielab::domain::se & a, const Lielab::domain::se & b, const int p)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
    
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

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("ad: Shapes of a and b must be equal.");
    }

    if (p == 0)
    {
        return b;
    }

    const ptrdiff_t shape = a.shape;
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
}

#endif
