#ifndef LIELAB_FUNCTIONS_EXTRA_H
#define LIELAB_FUNCTIONS_EXTRA_H

namespace Lielab::functions
{

// template<Lielab::abstract::LieAlgebra LA>
// Lielab::domain::cola coad(const LA & g, const Lielab::domain::cola & mu)
// {
//     /*! \f{equation}{ (\mathfrak{g}, \mathfrak{g}^*) \rightarrow \mathfrak{g}^* }
//     *
//     * Coadjoint function.
//     */

//     std::vector<LA> basis;
//     const size_t d1 = g.get_dimension();
//     const size_t d2 = mu.get_dimension();
//     // double norm = 0.0;

//     if (d1 != d2)
//     {
//         throw SizeError("Input data malformed: dimension " + std::to_string(d1) + " != " + std::to_string(d2) + ".");
//     }
    
//     Eigen::VectorXd out(d1);

//     const size_t s = g.shape;

//     for (int ii = 0; ii < d1; ii++)
//     {
//         basis.push_back(LA::basis(ii, s));
//     }

//     for (int ii = 0; ii < d1; ii++)
//     {
//         out(ii) = copair(mu, -ad(g, basis[ii]));
//     }

//     return out;
// }

}

#endif
