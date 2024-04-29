#ifndef _LIELAB_FUNCTIONS_CAYLEY1_HPP
#define _LIELAB_FUNCTIONS_CAYLEY1_HPP

namespace Lielab
{
    namespace functions
    {
        template <Lielab::abstract::LieAlgebra LA>
        Lielab::domain::lieiii<LA> cayley1(const LA & g)
        {
            /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
            *
            * Cayley transform.
            * 
            * \f{equation*}{ cayley1(g) = \frac{Id_g + g/2}{Id_g - g/2} \f}
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

            const Eigen::MatrixXd m = g.get_matrix();
            const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(g.shape, g.shape);

            return (Id + m/2.0)*(Id - m/2.0).inverse();
        }

        template <>
        Lielab::domain::CN cayley1(const Lielab::domain::cn & a)
        {
            /*
            * Cayley1 overload for cn
            *
            * Needed since cn is complex
            */

            const Eigen::MatrixXcd m = a.get_matrix();
            const Eigen::MatrixXcd Id = Eigen::MatrixXcd::Identity(a.shape, a.shape);

            return (Id + m/2.0)*(Id - m/2.0).inverse();
        }

        template <>
        Lielab::domain::GLC cayley1(const Lielab::domain::glc & a)
        {
            /*
            * Cayley1 overload for glc
            *
            * Needed since glc is complex
            */

            const Eigen::MatrixXcd m = a.get_matrix();
            const Eigen::MatrixXcd Id = Eigen::MatrixXcd::Identity(a.shape, a.shape);

            return (Id + m/2.0)*(Id - m/2.0).inverse();
        }

        template <>
        Lielab::domain::SU cayley1(const Lielab::domain::su & a)
        {
            /*
            * Cayley1 overload for su
            *
            * Needed since su is complex
            */

            const Eigen::MatrixXcd m = a.get_matrix();
            const Eigen::MatrixXcd Id = Eigen::MatrixXcd::Identity(a.shape, a.shape);

            return (Id + m/2.0)*(Id - m/2.0).inverse();
        }
    }
}

#endif
