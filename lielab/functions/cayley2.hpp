#ifndef _LIELAB_FUNCTIONS_CAYLEY2_HPP
#define _LIELAB_FUNCTIONS_CAYLEY2_HPP

namespace Lielab
{
    namespace functions
    {
        template <Lielab::abstract::LieAlgebra LA>
        Lielab::domain::lieiii<LA> cayley2(const LA & g)
        {
            /*! \f{equation*}{ (\mathfrak{g}) \rightarrow G \f}
            *
            * Cayley transform.
            * 
            * \f{equation*}{ cayley2(g) = cayley1(\xi_1 x_1)cayley1(\xi_2 x_2) \cdots cayley1(\xi_n x_n) \f}
            * 
            * @param[in] g A Lie algebra.
            * @param[out] G A Lie group.
            * 
            * Source: Engø, Kenth. "On the construction of geometric integrators in the RKMK class."
            * BIT Numerical Mathematics 40.1 (2000): 41-61.
            * 
            * TODO: Check that this works with SU. Math says no but simulations say otherwise.
            */

            const size_t dim = g.get_dimension();
            const size_t shape = g.shape;
            const Eigen::VectorXd v = g.get_vector();

            Lielab::domain::lieiii<LA> out(shape);

            for (int ii = 0; ii < dim; ii++)
            {
                const LA h = LA::basis(ii, shape);
                out = out*cayley1(v[ii]*h);
            }

            return out;
        }
    }
}

#endif
