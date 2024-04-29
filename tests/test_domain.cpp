#include <Lielab.hpp>
#include <iostream>

#include "test_utils.hpp"

using Lielab::functions::commutator;

/*!
* Checks each algebra up to the shape specified. Can get very expensive at
* large numbers
*/
constexpr size_t TEST_UP_TO_THIS_SHAPE = 4;


template <Lielab::abstract::Algebra T>
void is_algebra(const std::vector<T> & basis)
{
    /*!
    * Asserts whether or not a given basis forms an algebra.
    */
    
    const double a = 2.0;
    const double b = 3.0;

    for (auto & x : basis)
    {
        for (auto & y : basis)
        {
            // Scalar multiplication
            assert_domain<T>((a * x) * (b * y), (a * b) * (x * y));

            // Scalar division
            assert_domain<T>((x / a) * (y / b), (1.0 / (a * b)) * (x * y));

            // Vector addition
            assert_domain<T>(x + y, y + x);

            // Vector subtraction
            assert_domain<T>(x - y, -(y - x));

            for (auto & z : basis)
            {
                // Right distributive
                assert_domain<T>((x + y) * z, x * z + y * z);

                // Left distributive
                assert_domain<T>(x * (y + z), x * y + x * z);
            }
        }
    }
}

template <Lielab::abstract::Algebra T>
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
        // Alternating
        assert_domain<T>(commutator(x, x), zero);

        for (auto & y : basis)
        {
            // Anticommutivity
            assert_domain<T>(commutator(x, y), -commutator(y, x));

            // Abelian check
            if (x.abelian)
            {
                // Don't use commutator() since it will shortcut by returning 0.
                assert_domain<T>(x*y, y*x);
            }

            for (auto & z : basis)
            {
                // Bilinearity
                assert_domain<T>(commutator(a * x + b * y, z), a * commutator(x, z) + b * commutator(y, z));
                assert_domain<T>(commutator(z, a * x + b * y), a * commutator(z, x) + b * commutator(z, y));

                // Jacobi Identity
                assert_domain<T>(commutator(x, commutator(y, z)) + commutator(y, commutator(z, x)) + commutator(z, commutator(x, y)), zero);
            }
        }
    }
}

template <Lielab::abstract::Group T>
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


TEST_CASE("cn algebra", "[domain]")
{
    /*!
    * Tests the cn algebra.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::cn::basis(0, shape).get_dimension();

        // Construct the rn basis
        std::vector<Lielab::domain::cn> basis;
        for (size_t ii = 0; ii < D; ii++)
        {
            basis.push_back(Lielab::domain::cn::basis(ii, shape));
        }

        is_algebra<Lielab::domain::cn>(basis);
        is_liealgebra<Lielab::domain::cn>(basis);
    }
   
    Lielab::domain::cn one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(one.get_dimension() == 0);
    CHECK(two.get_dimension() == 2);
    CHECK(three.get_dimension() == 4);
    CHECK(four.get_dimension() == 6);
    CHECK(five.get_dimension() == 8);
    CHECK(six.get_dimension() == 10);
    CHECK(seven.get_dimension() == 12);
    CHECK(eight.get_dimension() == 14);
}


TEST_CASE("gl algebra", "[domain]")
{
    /*!
    * Tests the gl algebra.
    */

    for (size_t shape = 1; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::gl::basis(0, shape).get_dimension();

        // Construct the gl basis
        std::vector<Lielab::domain::gl> basis;
        for (size_t ii = 0; ii < D; ii++)
        {
            basis.push_back(Lielab::domain::gl::basis(ii, shape));
        }

        is_algebra<Lielab::domain::gl>(basis);
        is_liealgebra<Lielab::domain::gl>(basis);
    }
   
    Lielab::domain::gl one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(one.get_dimension() == 1);
    CHECK(two.get_dimension() == 4);
    CHECK(three.get_dimension() == 9);
    CHECK(four.get_dimension() == 16);
    CHECK(five.get_dimension() == 25);
    CHECK(six.get_dimension() == 36);
    CHECK(seven.get_dimension() == 49);
    CHECK(eight.get_dimension() == 64);
}


TEST_CASE("glc algebra", "[domain]")
{
    /*!
    * Tests the glc algebra.
    */

    for (size_t shape = 1; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::glc::basis(0, shape).get_dimension();

        // Construct the glc basis
        std::vector<Lielab::domain::glc> basis;
        for (size_t ii = 0; ii < D; ii++)
        {
            basis.push_back(Lielab::domain::glc::basis(ii, shape));
        }

        is_algebra<Lielab::domain::glc>(basis);
        is_liealgebra<Lielab::domain::glc>(basis);
    }
   
    Lielab::domain::glc one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(one.get_dimension() == 2);
    CHECK(two.get_dimension() == 8);
    CHECK(three.get_dimension() == 18);
    CHECK(four.get_dimension() == 32);
    CHECK(five.get_dimension() == 50);
    CHECK(six.get_dimension() == 72);
    CHECK(seven.get_dimension() == 98);
    CHECK(eight.get_dimension() == 128);
}


TEST_CASE("rn algebra", "[domain]")
{
    /*!
    * Tests the rn algebra.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::rn::basis(0, shape).get_dimension();

        // Construct the rn basis
        std::vector<Lielab::domain::rn> basis;
        for (size_t ii = 0; ii < D; ii++)
        {
            basis.push_back(Lielab::domain::rn::basis(ii, shape));
        }

        is_algebra<Lielab::domain::rn>(basis);
        is_liealgebra<Lielab::domain::rn>(basis);
    }
   
    Lielab::domain::rn one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(one.get_dimension() == 0);
    CHECK(two.get_dimension() == 1);
    CHECK(three.get_dimension() == 2);
    CHECK(four.get_dimension() == 3);
    CHECK(five.get_dimension() == 4);
    CHECK(six.get_dimension() == 5);
    CHECK(seven.get_dimension() == 6);
    CHECK(eight.get_dimension() == 7);
}

TEST_CASE("se algebra", "[domain]")
{
    /*!
    * Tests the se algebra.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::se::basis(0, shape).get_dimension();

        // Construct the so basis
        std::vector<Lielab::domain::se> basis;
        for (size_t ii = 0; ii < D; ii++)
        {
            basis.push_back(Lielab::domain::se::basis(ii, shape));
        }

        is_algebra<Lielab::domain::se>(basis);
        is_liealgebra<Lielab::domain::se>(basis);
    }
   
    Lielab::domain::se one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(one.get_dimension() == 0);
    CHECK(two.get_dimension() == 1);
    CHECK(three.get_dimension() == 3);
    CHECK(four.get_dimension() == 6);
    CHECK(five.get_dimension() == 10);
    CHECK(six.get_dimension() == 15);
    CHECK(seven.get_dimension() == 21);
    CHECK(eight.get_dimension() == 28);
}

TEST_CASE("se4", "[domain]")
{
    /*!
    * Tests the se algebra with se(4).
    */

    Lielab::domain::se x = Lielab::domain::se::basis(0, 4);
    Lielab::domain::se y = Lielab::domain::se::basis(1, 4);
    Lielab::domain::se z = Lielab::domain::se::basis(2, 4);
    Lielab::domain::se u = Lielab::domain::se::basis(3, 4);
    Lielab::domain::se v = Lielab::domain::se::basis(4, 4);
    Lielab::domain::se w = Lielab::domain::se::basis(5, 4);
    Lielab::domain::se zero = x*0;

    assert_domain(commutator(x, y), zero);
    assert_domain(commutator(y, z), zero);
    assert_domain(commutator(z, x), zero);
    assert_domain(commutator(y, x), zero);
    assert_domain(commutator(z, y), zero);
    assert_domain(commutator(x, z), zero);
    assert_domain(commutator(x, v), z);
    assert_domain(commutator(y, w), x);
    assert_domain(commutator(z, u), y);
    assert_domain(commutator(y, u), -z);
    assert_domain(commutator(z, v), -x);
    assert_domain(commutator(x, w), -y);
    assert_domain(commutator(u, v), w);
    assert_domain(commutator(v, w), u);
    assert_domain(commutator(w, u), v);
    assert_domain(commutator(v, u), -w);
    assert_domain(commutator(w, v), -u);
    assert_domain(commutator(u, w), -v);
}

TEST_CASE("so algebra", "[domain]")
{
    /*!
    * Tests the so algebra.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::so::basis(0, shape).get_dimension();

        // Construct the so basis
        std::vector<Lielab::domain::so> basis;
        for (size_t ii = 0; ii < D; ii++)
        {
            basis.push_back(Lielab::domain::so::basis(ii, shape));
        }

        is_algebra<Lielab::domain::so>(basis);
        is_liealgebra<Lielab::domain::so>(basis);
    }
   
    Lielab::domain::so one(1), two(2), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(one.get_dimension() == 0);
    CHECK(two.get_dimension() == 1);
    // so3 is checked independently
    CHECK(four.get_dimension() == 6);
    CHECK(five.get_dimension() == 10);
    CHECK(six.get_dimension() == 15);
    CHECK(seven.get_dimension() == 21);
    CHECK(eight.get_dimension() == 28);
}

TEST_CASE("so3", "[domain]")
{
    /*!
    * Tests the so algebra with so(3).
    */

    Lielab::domain::so x = Lielab::domain::so::basis(0,3);
    Lielab::domain::so y = Lielab::domain::so::basis(1,3);
    Lielab::domain::so z = Lielab::domain::so::basis(2,3);
    Lielab::domain::so zero = x*0;

    assert_domain(commutator(x, y), z);
    assert_domain(commutator(y, z), x);
    assert_domain(commutator(z, x), y);
    assert_domain(commutator(y, x), -z);
    assert_domain(commutator(z, y), -x);
    assert_domain(commutator(x, z), -y);
}

TEST_CASE("sp algebra", "[domain]")
{
    /*!
    * Tests the sp algebra.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape += 2)
    {
        const size_t D = Lielab::domain::sp::basis(0, shape).get_dimension();

        // Construct the sp basis
        std::vector<Lielab::domain::sp> basis;
        for (size_t ii = 0; ii < D; ii++)
        {
            basis.push_back(Lielab::domain::sp::basis(ii, shape));
        }

        is_algebra<Lielab::domain::sp>(basis);
        is_liealgebra<Lielab::domain::sp>(basis);
    }
   
    Lielab::domain::sp two(2), four(4), six(6), eight(8);

    // Dimensions
    CHECK(two.get_dimension() == 3);
    CHECK(four.get_dimension() == 10);
    CHECK(six.get_dimension() == 21);
    CHECK(eight.get_dimension() == 36);
}

TEST_CASE("su", "[domain]")
{
    /*!
    * Tests the su algebra implementation.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::su::basis(0, shape).get_dimension();

        // Construct the su basis
        std::vector<Lielab::domain::su> basis;
        for (size_t ii = 0; ii < D; ii++)
        {
            basis.push_back(Lielab::domain::su::basis(ii, shape));
        }

        is_algebra<Lielab::domain::su>(basis);
        is_liealgebra<Lielab::domain::su>(basis);
    }
}

TEST_CASE("su2", "[domain]")
{
    /*!
    * Tests the su algebra with su(2).
    */

    const size_t D = Lielab::domain::su::basis(0, 2).get_dimension();

    // Construct the su2 basis
    std::vector<Lielab::domain::su> b;
    for (size_t ii = 0; ii < D; ii++)
    {
        b.push_back(Lielab::domain::su::basis(ii, 2));
    }

    // su2 specific identities
    assert_domain(commutator(b[0], b[1]),  2*b[2]);
    assert_domain(commutator(b[1], b[2]),  2*b[0]);
    assert_domain(commutator(b[2], b[0]),  2*b[1]);
    assert_domain(commutator(b[1], b[0]), -2*b[2]);
    assert_domain(commutator(b[2], b[1]), -2*b[0]);
    assert_domain(commutator(b[0], b[2]), -2*b[1]);

    // Hamilton's identities
    // Note that i^2 = j^2 = k^2 = -1^2 isn't checked since this isn't true for the algebra

    // ij = -ji = k
    assert_domain(b[0]*b[1], -b[1]*b[0]);
    assert_domain(b[0]*b[1],  b[2]);

    // jk = -kj = i
    assert_domain(b[1]*b[2], -b[2]*b[1]);
    assert_domain(b[1]*b[2],  b[0]);

    // ki = -ik = j
    assert_domain(b[2]*b[0], -b[0]*b[2]);
    assert_domain(b[2]*b[0],  b[1]);
}

TEST_CASE("su3", "[domain]")
{
    /*!
    * Tests the su algebra with su(3).
    */

    const size_t D = Lielab::domain::su::basis(0, 3).get_dimension();

    // Construct the su3 basis
    std::vector<Lielab::domain::su> b;
    for (size_t ii = 0; ii < D; ii++)
    {
        b.push_back(Lielab::domain::su::basis(ii, 3));
    }

    // TODO: Implement these. It's tough to work out what these should be with GGM
    // su3 specific identities
    // assert_domain( commutator(t1, t2), _i*t3);
    // assert_domain( commutator(t1, t4), _i*t7 / 2.0);
    // assert_domain(-commutator(t1, t5), _i*t6 / 2.0);
    // assert_domain( commutator(t2, t4), _i*t6 / 2.0);
    // assert_domain( commutator(t2, t5), _i*t7 / 2.0);
    // assert_domain( commutator(t3, t4), _i*t5 / 2.0);
    // assert_domain(-commutator(t3, t6), _i*t7 / 2.0);
    // TODO: Why aren't these evaluating true?
    // assert_domain( commutator(t4, t5), _i*t8 * std::sqrt(3.0) / 2.0);
    // assert_domain( commutator(t6, t7), _i*t8 * std::sqrt(3.0) / 2.0);
}

TEST_CASE("CN", "[domain]")
{
    /*!
    * Tests CN against well-known identities.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::cn::basis(0, shape).get_dimension();

        // Construct the RN elements
        std::vector<Lielab::domain::CN> elements;
        for (size_t ii = 0; ii < D; ii++)
        {
            elements.push_back(Lielab::functions::exp(Lielab::domain::cn::basis(ii, shape)));
        }

        const Lielab::domain::CN identity(shape);

        is_group<Lielab::domain::CN>(elements, identity);
    }
}

TEST_CASE("GL", "[domain]")
{
    /*!
    * Tests GL against well-known identities.
    */

    for (size_t shape = 1; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::gl::basis(0, shape).get_dimension();

        // Construct the GL elements
        std::vector<Lielab::domain::GL> elements;
        for (size_t ii = 0; ii < D; ii++)
        {
            elements.push_back(Lielab::functions::exp(Lielab::domain::gl::basis(ii, shape)));
        }

        const Lielab::domain::GL identity(shape);

        is_group<Lielab::domain::GL>(elements, identity);
    }
}

TEST_CASE("GLC", "[domain]")
{
    /*!
    * Tests GLC against well-known identities.
    */

    for (size_t shape = 1; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::glc::basis(0, shape).get_dimension();

        // Construct the GL elements
        std::vector<Lielab::domain::GLC> elements;
        for (size_t ii = 0; ii < D; ii++)
        {
            elements.push_back(Lielab::functions::exp(Lielab::domain::glc::basis(ii, shape)));
        }

        const Lielab::domain::GLC identity(shape);

        is_group<Lielab::domain::GLC>(elements, identity);
    }
}

TEST_CASE("RN", "[domain]")
{
    /*!
    * Tests RN against well-known identities.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::rn::basis(0, shape).get_dimension();

        // Construct the RN elements
        std::vector<Lielab::domain::RN> elements;
        for (size_t ii = 0; ii < D; ii++)
        {
            elements.push_back(Lielab::functions::exp(Lielab::domain::rn::basis(ii, shape)));
        }

        const Lielab::domain::RN identity(shape);

        is_group<Lielab::domain::RN>(elements, identity);
    }
}

TEST_CASE("SO", "[domain]")
{
    /*!
    * Tests SO against well-known identities.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::so::basis(0, shape).get_dimension();

        // Construct the SO elements
        std::vector<Lielab::domain::SO> elements;
        for (size_t ii = 0; ii < D; ii++)
        {
            elements.push_back(Lielab::functions::exp(Lielab::domain::so::basis(ii, shape)));
        }

        const Lielab::domain::SO identity(shape);

        is_group<Lielab::domain::SO>(elements, identity);
    }
}

TEST_CASE("SP", "[domain]")
{
    /*!
    * Tests SP against well-known identities.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape += 2)
    {
        const size_t D = Lielab::domain::sp::basis(0, shape).get_dimension();

        // Construct the SP elements
        std::vector<Lielab::domain::SP> elements;
        for (size_t ii = 0; ii < D; ii++)
        {
            elements.push_back(Lielab::functions::exp(Lielab::domain::sp::basis(ii, shape)));
        }

        const Lielab::domain::SP identity(shape);

        is_group<Lielab::domain::SP>(elements, identity);
    }
}

TEST_CASE("SU", "[domain]")
{
    /*!
    * Tests SU against well-known identities.
    */

    for (size_t shape = 2; shape <= TEST_UP_TO_THIS_SHAPE; shape++)
    {
        const size_t D = Lielab::domain::su::basis(0, shape).get_dimension();

        // Construct the SU elements
        std::vector<Lielab::domain::SU> elements;
        for (size_t ii = 0; ii < D; ii++)
        {
            elements.push_back(Lielab::functions::exp(Lielab::domain::su::basis(ii, shape)));
        }

        const Lielab::domain::SU identity(shape);

        is_group<Lielab::domain::SU>(elements, identity);
    }
}

TEST_CASE("Quaternion", "[domain]")
{
    /*!
    * Tests quaternions against well-known identities.
    */

    Lielab::domain::SU q1 = Lielab::domain::SU::Quaternion();
    Lielab::domain::SU qm1 = Lielab::domain::SU::Quaternion(-1.0, 0.0, 0.0, 0.0);
    Lielab::domain::SU qi = Lielab::domain::SU::Quaternion(0.0, 1.0, 0.0, 0.0);
    Lielab::domain::SU qj = Lielab::domain::SU::Quaternion(0.0, 0.0, 1.0, 0.0);
    Lielab::domain::SU qk = Lielab::domain::SU::Quaternion(0.0, 0.0, 0.0, 1.0);

    std::vector<Lielab::domain::SU> elements;
    elements.push_back(qi);
    elements.push_back(qj);
    elements.push_back(qk);

    is_group<Lielab::domain::SU>(elements, q1);

    // Hamilton's identities
    // i^2 = j^2 = k^2 = -1
    assert_domain(qi*qi, qm1);
    assert_domain(qj*qj, qm1);
    assert_domain(qk*qk, qm1);

    // ij = -ji = -k
    assert_domain(qi*qj, (qj*qi).inverse());
    assert_domain(qi*qj, qk.inverse());

    // jk = -kj = -i
    assert_domain(qj*qk, (qk*qj).inverse());
    assert_domain(qj*qk, qi.inverse());

    // ki = -ik = -j
    assert_domain(qk*qi, (qi*qk).inverse());
    assert_domain(qk*qi, qj.inverse());
}
