#ifndef _LIELAB_FUNCTIONS_KILLING_HPP
#define _LIELAB_FUNCTIONS_KILLING_HPP

namespace Lielab
{
    namespace functions
    {
        template <Lielab::abstract::LieAlgebra LA>
        double Killing(LA & g, LA & h)
        {
            /*!
            * 
            * @param[in] g Lie algebra g.
            * @param[in] h Lie algebra h.
            * @param[out] k Killing coefficient between g and h.
            */

            const size_t dim = g.get_dimension();
            const size_t shape = g.shape;

            std::vector<LA> basis = std::vector<LA>(dim);

            for (int ii = 0; ii < dim; ii++)
            {
                basis[ii] = LA::basis(ii, shape);
            }

            double k = 0.0;
            for (int ii = 0; ii < dim; ii++)
            {
                k += commutator(g, commutator(h, basis[ii])).get_vector()[ii];
            }

            return k;
        }
    }
}

#endif
