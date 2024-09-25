#ifndef _LIELAB_FUNCTIONS_KILLING_HPP
#define _LIELAB_FUNCTIONS_KILLING_HPP

#include "../abstract.hpp"
#include "../domain.hpp"

namespace Lielab
{
namespace functions
{
template <Lielab::abstract::LieAlgebra LA>
double Killing(const LA & a, const LA & b)
{
    /*!
    * 
    * @param[in] a First Lie algebra element.
    * @param[in] b Second Lie algebra element
    * @param[out] k Killing coefficient between a and b.
    */

    if (a.shape != b.shape)
    {
        throw Lielab::utils::InputError("Killing: Shapes of a and b must be equal.");
    }

    const size_t dim = a.get_dimension();
    const size_t shape = a.shape;

    std::vector<LA> basis = std::vector<LA>(dim);

    for (int ii = 0; ii < dim; ii++)
    {
        basis[ii] = LA::basis(ii, shape);
    }

    double k = 0.0;
    for (int ii = 0; ii < dim; ii++)
    {
        k += commutator(a, commutator(b, basis[ii])).get_vector()[ii];
    }

    return k;
}
}
}

#endif
