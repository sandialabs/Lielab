#ifndef LIELAB_FUNCTIONS_KILLING_TPP
#define LIELAB_FUNCTIONS_KILLING_TPP

#include "Killing.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <typename LA>
double Killing(const LA& a, const LA& b)
{
    /*!
    * 
    * @param[in] a First Lie algebra element.
    * @param[in] b Second Lie algebra element
    * @param[out] k Killing coefficient between a and b.
    */

    const int shape = a.get_shape();

    if (shape != b.get_shape())
    {
        throw Lielab::utils::InputError("Killing: Shapes of a and b must be equal.");
    }

    const int dim = a.get_dimension();

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

template <typename LA>
Eigen::MatrixXd Killingform(const LA& g)
{
    /*!
    *
    * @param[in] g Lie algebra g.
    * @param[out] K Killing form of g in quadratic / matrix form.
    * TODO: Make this a class method?
    * TODO: Return as glr type?
    */

    const int dim = g.get_dimension();
    const int shape = g.get_shape();
    std::vector<LA> basis = std::vector<LA>(dim);

    for (int ii = 0; ii < dim; ii++)
    {
        basis[ii] = LA::basis(ii, shape);
    }

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(dim, dim);

    for (int ii = 0; ii < dim; ii++)
    {
        for (int jj = 0; jj < dim; jj++)
        {
            K(ii, jj) = Killing(basis[ii], basis[jj]);
        }
    }

    return K;
}

}

#endif
