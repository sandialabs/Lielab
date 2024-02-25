#ifndef _LIELAB_DOMAIN_KILLINGFORM_HPP
#define _LIELAB_DOMAIN_KILLINGFORM_HPP

namespace lielab
{
    namespace functions
    {
        template <lielab::abstract::LieAlgebra LA>
        Eigen::MatrixXd Killingform(LA & g)
        {
            /*!
            *
            * @param[in] g Lie algebra g.
            * @param[out] K Killing form of g in quadratic / matrix form.
            * TODO: Make this a class method?
            */

            const size_t dim = g.get_dimension();
            const size_t shape = g.shape;
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
}

#endif
