#ifndef LIELAB_TEST_TPP
#define LIELAB_TEST_TPP

#include <iostream>

#include <Eigen/Core>
#include <Lielab.hpp>

template<typename L>
void assert_domain(L mat1, L mat2)
{
    assert_matrix(mat1.get_matrix(), mat2.get_matrix());
}

template <typename T>
void is_algebra(const std::vector<T> & basis)
{
    /*!
    * Asserts whether or not a given basis forms an algebra.
    */
    
    const double a = 2.0;
    const double b = 3.0;

    for (auto & x : basis)
    {
        const auto xhat = x.get_matrix();

        for (auto & y : basis)
        {
            const auto yhat = y.get_matrix();

            // Scalar multiplication
            assert_domain<T>((a * xhat) * (b * yhat), (a * b) * (xhat * yhat));

            // Scalar division
            assert_domain<T>((xhat / a) * (yhat / b), (1.0 / (a * b)) * (xhat * yhat));

            // Vector addition
            assert_domain<T>(x + y, y + x);

            // Vector subtraction
            assert_domain<T>(x - y, -(y - x));

            for (auto & z : basis)
            {
                const auto zhat = z.get_matrix();

                // Right distributive
                assert_domain<T>((xhat + yhat) * zhat, xhat * zhat + yhat * zhat);

                // Left distributive
                assert_domain<T>(xhat * (yhat + zhat), xhat * yhat + xhat * zhat);
            }
        }
    }
}

template <typename T>
void is_liealgebra(const std::vector<T> & basis)
{
    /*!
    * Asserts whether or not a given basis forms a Lie algebra.
    */

    const double a = 2.0;
    const double b = 3.0;
    
    const T zero = (basis.size() > 0) ? 0.0*basis[0] : T(0);

    for (auto & x : basis)
    {
        const auto xhat = x.get_matrix();

        // Alternating
        assert_domain<T>(Lielab::functions::commutator(x, x), zero);

        for (auto & y : basis)
        {
            const auto yhat = y.get_matrix();

            // Anticommutivity
            assert_domain<T>(Lielab::functions::commutator(x, y), -Lielab::functions::commutator(y, x));

            // Abelian check
            if (x.abelian)
            {
                // Don't use commutator() since it will shortcut by returning 0.
                assert_domain<T>(xhat*yhat, yhat*xhat);
            }

            for (auto & z : basis)
            {
                // Bilinearity
                assert_domain<T>(Lielab::functions::commutator(a * x + b * y, z), a * Lielab::functions::commutator(x, z) + b * Lielab::functions::commutator(y, z));
                assert_domain<T>(Lielab::functions::commutator(z, a * x + b * y), a * Lielab::functions::commutator(z, x) + b * Lielab::functions::commutator(z, y));

                // Jacobi Identity
                assert_domain<T>(Lielab::functions::commutator(x, Lielab::functions::commutator(y, z)) + Lielab::functions::commutator(y, Lielab::functions::commutator(z, x)) + Lielab::functions::commutator(z, Lielab::functions::commutator(x, y)), zero);
            }
        }
    }
}

template <typename T>
void is_group(const std::vector<T> & elements, const T & identity)
{
    /*!
    * Asserts whether or not a given set of elements are in a group.
    */

    for (auto & x : elements)
    {
        // Identity
        assert_domain<T>(x * identity, x);
        assert_domain<T>(identity * x, x);

        // Inverse
        assert_domain<T>(x * x.inverse(), identity);

        for (auto & y : elements)
        {
            // Inverse
            assert_domain<T>(x*y, (y.inverse() * x.inverse()).inverse());

            // Abelian check
            if (x.abelian)
            {
                assert_domain<T>(x*y, y*x);
            }
            
            for (auto & z : elements)
            {
                // Associative
                assert_domain<T>((x * y) * z, x * (y * z));
            }
        }
    }
}

#endif
