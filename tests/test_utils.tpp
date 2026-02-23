#ifndef LIELAB_TEST_TPP
#define LIELAB_TEST_TPP

#include <iostream>

#include <Eigen/Core>
#include <Lielab.hpp>

template <typename T>
void is_liealgebra(const std::vector<T>& basis)
{
    /*!
    * Asserts whether or not a given basis forms a Lie algebra.
    */

    using Lielab::functions::commutator;
    using Lielab::testing::check_almost_equal_nulp;

    const double a = 2.0;
    const double b = 3.0;
    
    const T zero = (basis.size() > 0) ? 0.0*basis[0] : T(0);

    for (auto& x : basis)
    {
        const auto xhat = x.get_matrix();

        // Alternating
        CHECK(check_almost_equal_nulp(commutator(x, x).get_matrix(), zero.get_matrix(), 1, true));

        for (auto& y : basis)
        {
            const auto yhat = y.get_matrix();

            // Anticommutivity
            CHECK(check_almost_equal_nulp(commutator(x, y).get_matrix(), (-commutator(y, x)).get_matrix(), 1, true));

            // Abelian check
            if (x.abelian)
            {
                // Don't use commutator() since it will shortcut by returning 0.
                CHECK(check_almost_equal_nulp((xhat*yhat).get_matrix(), (yhat*xhat).get_matrix(), 1, true));
            }

            for (auto& z : basis)
            {
                // Bilinearity
                CHECK(check_almost_equal_nulp(commutator(a * x + b * y, z).get_matrix(), (a * commutator(x, z) + b * commutator(y, z)).get_matrix(), 1, true));
                CHECK(check_almost_equal_nulp(commutator(z, a * x + b * y).get_matrix(), (a * commutator(z, x) + b * commutator(z, y)).get_matrix(), 1, true));

                // Jacobi Identity
                CHECK(check_almost_equal_nulp((commutator(x, commutator(y, z)) + commutator(y, commutator(z, x)) + commutator(z, commutator(x, y))).get_matrix(), zero.get_matrix(), 1, true));
            }
        }
    }
}

template <typename T>
void is_group(const std::vector<T>& elements, const T& identity)
{
    /*!
    * Asserts whether or not a given set of elements are in a group.
    */

    using Lielab::testing::check_almost_equal_nulp;

    for (auto& x : elements)
    {
        // Identity
        CHECK(check_almost_equal_nulp((x * identity).get_matrix(), x.get_matrix(), 1, true));
        CHECK(check_almost_equal_nulp((identity * x).get_matrix(), x.get_matrix(), 1, true));

        // Inverse
        CHECK(check_almost_equal_nulp((x * x.inverse()).get_matrix(), identity.get_matrix(), 1, true));

        for (auto& y : elements)
        {
            // Inverse
            CHECK(check_almost_equal_nulp((x*y).get_matrix(), (y.inverse() * x.inverse()).inverse().get_matrix(), 1, true));

            // Abelian check
            if (x.is_abelian())
            {
                CHECK(check_almost_equal_nulp((x*y).get_matrix(), (y*x).get_matrix(), 1, true));
            }
            
            for (auto& z : elements)
            {
                // Associative
                CHECK(check_almost_equal_nulp(((x * y) * z).get_matrix(), (x * (y * z)).get_matrix(), 1, true));
            }
        }
    }
}

#endif
