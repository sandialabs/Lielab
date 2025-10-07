#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

#include <iostream>

#include "../../test_utils.hpp"

TEST_CASE("su to_string", "[domain]")
{
    using namespace Lielab::domain;

    const su x0 = su(0);
    CHECK(x0.to_string() == "su(0)");
    const su x1 = su(1);
    CHECK(x1.to_string() == "su(1)");
    const su x10 = su(10);
    CHECK(x10.to_string() == "su(10)");
}

TEST_CASE("su main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const su xblank = su();
    CHECK(xblank.get_dimension() == 0);

    const su x0 = su(0);
    CHECK(x0.get_dimension() == 0);
    const su x1 = su(1);
    CHECK(x1.get_dimension() == 0);
    const su x10 = su(2);
    CHECK(x10.get_dimension() == 3);
}

TEST_CASE("su matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const su x0 = su(Eigen::MatrixXcd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    const su x1 = su(Eigen::MatrixXcd::Random(1, 1));
    CHECK(x1.get_shape() == 1);

    const su x2 = su(Eigen::MatrixXcd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(su(Eigen::MatrixXcd::Random(2, 3)));
    CHECK_THROWS(su(Eigen::MatrixXcd::Random(3, 2)));
}

TEST_CASE("su basis_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const su xm10 = su::basis(-1, 0);
    CHECK(xm10.get_dimension() == 0);
    const Eigen::VectorXd xm10bar = xm10.get_vector();
    CHECK(xm10bar.size() == 0);

    const su x00 = su::basis(0, 0);
    CHECK(x00.get_dimension() == 0);
    const Eigen::VectorXd x00bar = x00.get_vector();
    CHECK(x00bar.size() == 0);

    const su x10 = su::basis(1, 0);
    CHECK(x10.get_dimension() == 0);
    const Eigen::VectorXd x10bar = x10.get_vector();
    CHECK(x10bar.size() == 0);

    const su x02 = su::basis(0, 2);
    CHECK(x02.get_dimension() == 3);
    const Eigen::VectorXd x02bar = x02.get_vector();
    REQUIRE(x02bar.size() == 3);
    CHECK(x02bar(0) == 1.0);
    CHECK(x02bar(1) == 0.0);
    CHECK(x02bar(2) == 0.0);

    const su x12 = su::basis(1, 2);
    CHECK(x12.get_dimension() == 3);
    const Eigen::VectorXd x12bar = x12.get_vector();
    REQUIRE(x12bar.size() == 3);
    CHECK(x12bar(0) == 0.0);
    CHECK(x12bar(1) == 1.0);
    CHECK(x12bar(2) == 0.0);

    const su x22 = su::basis(2, 2);
    CHECK(x22.get_dimension() == 3);
    const Eigen::VectorXd x22bar = x22.get_vector();
    REQUIRE(x22bar.size() == 3);
    CHECK(x22bar(0) == 0.0);
    CHECK(x22bar(1) == 0.0);
    CHECK(x22bar(2) == 1.0);

    const su x32 = su::basis(3, 2);
    CHECK(x32.get_dimension() == 3);
    const Eigen::VectorXd x32bar = x32.get_vector();
    REQUIRE(x32bar.size() == 3);
    CHECK(x32bar(0) == 0.0);
    CHECK(x32bar(1) == 0.0);
    CHECK(x32bar(2) == 0.0);

    const su x03 = su::basis(0, 3);
    CHECK(x03.get_dimension() == 8);
    const Eigen::VectorXd x03bar = x03.get_vector();
    REQUIRE(x03bar.size() == 8);
    CHECK(x03bar(0) == 1.0);
    CHECK(x03bar(1) == 0.0);
    CHECK(x03bar(2) == 0.0);
    CHECK(x03bar(3) == 0.0);
    CHECK(x03bar(4) == 0.0);
    CHECK(x03bar(5) == 0.0);
    CHECK(x03bar(6) == 0.0);
    CHECK(x03bar(7) == 0.0);
}

TEST_CASE("su from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const su x0 = su::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXcd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const su x1 = su::from_shape(1);
    CHECK(x1.get_dimension() == 0);
    const Eigen::MatrixXcd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const su x2 = su::from_shape(2);
    CHECK(x2.get_dimension() == 3);
    const Eigen::MatrixXcd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("su get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    su zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    CHECK(zero.get_dimension() == 0);
    CHECK(one.get_dimension() == 0);
    CHECK(two.get_dimension() == 3);
    CHECK(three.get_dimension() == 8);
    CHECK(four.get_dimension() == 15);
    CHECK(five.get_dimension() == 24);
    CHECK(six.get_dimension() == 35);
    CHECK(seven.get_dimension() == 48);
    CHECK(eight.get_dimension() == 63);
}

TEST_CASE("su set/get_vector", "[domain]")
{
    /*!
    * Tests the set/get_vector operation.
    */

    using namespace Lielab::domain;

    su x0 = su(0);
    x0.set_vector({});
    Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0bar.size() == 0);

    su x2 = su(2);
    x2.set_vector({1.0, 2.0, 3.0, 4.0});
    Eigen::VectorXd x2bar = x2.get_vector();

    REQUIRE(x2bar.size() == 3);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 3.0);

    x2.set_vector({5.0, 6.0, 7.0});
    x2bar = x2.get_vector();

    REQUIRE(x2bar.size() == 3);
    CHECK(x2bar(0) == 5.0);
    CHECK(x2bar(1) == 6.0);
    CHECK(x2bar(2) == 7.0);

    x2.set_vector({8.0, 9.0});
    x2bar = x2.get_vector();

    REQUIRE(x2bar.size() == 3);
    CHECK(x2bar(0) == 8.0);
    CHECK(x2bar(1) == 9.0);
    CHECK(x2bar(2) == 7.0);

    su x3 = su(3);
    x3.set_vector({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});
    Eigen::VectorXd x3bar = x3.get_vector();

    REQUIRE(x3bar.size() == 8);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
    CHECK(x3bar(3) == 4.0);
    CHECK(x3bar(4) == 5.0);
    CHECK(x3bar(5) == 6.0);
    CHECK(x3bar(6) == 7.0);
    CHECK(x3bar(7) == 0.0);

    x3.set_vector({8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0});
    x3bar = x3.get_vector();

    REQUIRE(x3bar.size() == 8);
    CHECK(x3bar(0) == 8.0);
    CHECK(x3bar(1) == 9.0);
    CHECK(x3bar(2) == 10.0);
    CHECK(x3bar(3) == 11.0);
    CHECK(x3bar(4) == 12.0);
    CHECK(x3bar(5) == 13.0);
    CHECK(x3bar(6) == 14.0);
    CHECK(x3bar(7) == 15.0);

    x3.set_vector({16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0});
    x3bar = x3.get_vector();

    REQUIRE(x3bar.size() == 8);
    CHECK(x3bar(0) == 16.0);
    CHECK(x3bar(1) == 17.0);
    CHECK(x3bar(2) == 18.0);
    CHECK(x3bar(3) == 19.0);
    CHECK(x3bar(4) == 20.0);
    CHECK(x3bar(5) == 21.0);
    CHECK(x3bar(6) == 22.0);
    CHECK(x3bar(7) == 23.0);
}

TEST_CASE("su get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    su x0 = su(0);
    x0.set_vector({});
    Eigen::MatrixXcd x0hat = x0.get_matrix();

    REQUIRE(x0hat.rows() == 0);
    REQUIRE(x0hat.cols() == 0);

    su x1 = su(1);
    x1.set_vector({});
    Eigen::MatrixXcd x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == std::complex<double>(0.0, 0.0));

    su x2 = su(2);
    x2.set_vector({1.0, 2.0, 3.0});
    Eigen::MatrixXcd x2hat = x2.get_matrix();

    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == std::complex<double>(0.0, 3.0));
    CHECK(x2hat(0, 1) == std::complex<double>(-2.0, 1.0));
    CHECK(x2hat(1, 0) == std::complex<double>(2.0, 1.0));
    CHECK(x2hat(1, 1) == std::complex<double>(0.0, -3.0));

    x2.set_vector({4.0, 5.0});
    x2hat = x2.get_matrix();

    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == std::complex<double>(0.0, 3.0));
    CHECK(x2hat(0, 1) == std::complex<double>(-5.0, 4.0));
    CHECK(x2hat(1, 0) == std::complex<double>(5.0, 4.0));
    CHECK(x2hat(1, 1) == std::complex<double>(0.0, -3.0));

    x2.set_vector({6.0});
    x2hat = x2.get_matrix();

    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == std::complex<double>(0.0, 3.0));
    CHECK(x2hat(0, 1) == std::complex<double>(-5.0, 6.0));
    CHECK(x2hat(1, 0) == std::complex<double>(5.0, 6.0));
    CHECK(x2hat(1, 1) == std::complex<double>(0.0, -3.0));
}

TEST_CASE("su operator()", "[domain]")
{
    using namespace Lielab::domain;

    su x0 = su(0);
    x0.set_vector({});

    // Out of bounds
    CHECK(std::isnan(x0(-1)));
    CHECK(std::isnan(x0(0)));
    CHECK(std::isnan(x0(1)));

    // Out of bounds
    CHECK(std::isnan(x0(0, -1).real()));
    CHECK(std::isnan(x0(0, -1).imag()));
    CHECK(std::isnan(x0(-1, 0).real()));
    CHECK(std::isnan(x0(-1, 0).imag()));
    CHECK(std::isnan(x0(-1, -1).real()));
    CHECK(std::isnan(x0(-1, -1).imag()));
    CHECK(std::isnan(x0(0, 0).real()));
    CHECK(std::isnan(x0(0, 0).imag()));
    CHECK(std::isnan(x0(0, 1).real()));
    CHECK(std::isnan(x0(0, 1).imag()));
    CHECK(std::isnan(x0(1, 0).real()));
    CHECK(std::isnan(x0(1, 0).imag()));
    CHECK(std::isnan(x0(1, 1).real()));
    CHECK(std::isnan(x0(1, 1).imag()));

    su x1 = su(1);
    x1.set_vector({});

    // Out of bounds
    CHECK(std::isnan(x1(-1)));
    CHECK(std::isnan(x1(0)));
    CHECK(std::isnan(x1(1)));

    // In bounds
    CHECK(x1(0, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x1(-1, -1) == std::complex<double>(0.0, 0.0));

    // Out of bounds
    CHECK(std::isnan(x1(0, -2).real()));
    CHECK(std::isnan(x1(0, -2).imag()));
    CHECK(std::isnan(x1(-2, 0).real()));
    CHECK(std::isnan(x1(-2, 0).imag()));
    CHECK(std::isnan(x1(-2, -2).real()));
    CHECK(std::isnan(x1(-2, -2).imag()));
    CHECK(std::isnan(x1(0, 1).real()));
    CHECK(std::isnan(x1(0, 1).imag()));
    CHECK(std::isnan(x1(1, 0).real()));
    CHECK(std::isnan(x1(1, 0).imag()));
    CHECK(std::isnan(x1(1, 1).real()));
    CHECK(std::isnan(x1(1, 1).imag()));

    su x2 = su(2);
    x2.set_vector({1.0, 2.0, 3.0});

    // In bounds
    CHECK(x2(0) == 1.0);
    CHECK(x2(1) == 2.0);
    CHECK(x2(2) == 3.0);
    CHECK(x2(-1) == 3.0);
    CHECK(x2(-2) == 2.0);
    CHECK(x2(-3) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x2(-4)));
    CHECK(std::isnan(x2(3)));

    // In bounds
    CHECK(x2(0, 0) == std::complex<double>(0.0, 3.0));
    CHECK(x2(0, 1) == std::complex<double>(-2.0, 1.0));
    CHECK(x2(1, 0) == std::complex<double>(2.0, 1.0));
    CHECK(x2(1, 1) == std::complex<double>(0.0, -3.0));
    CHECK(x2(-1, -1) == std::complex<double>(0.0, -3.0));
    CHECK(x2(-1, -2) == std::complex<double>(2.0, 1.0));
    CHECK(x2(-2, -1) == std::complex<double>(-2.0, 1.0));
    CHECK(x2(-2, -2) == std::complex<double>(0.0, 3.0));

    // Out of bounds
    CHECK(std::isnan(x2(0, -3).real()));
    CHECK(std::isnan(x2(0, -3).imag()));
    CHECK(std::isnan(x2(-3, 0).real()));
    CHECK(std::isnan(x2(-3, 0).imag()));
    CHECK(std::isnan(x2(-3, -3).real()));
    CHECK(std::isnan(x2(-3, -3).imag()));
    CHECK(std::isnan(x2(0, 2).real()));
    CHECK(std::isnan(x2(0, 2).imag()));
    CHECK(std::isnan(x2(2, 0).real()));
    CHECK(std::isnan(x2(2, 0).imag()));
    CHECK(std::isnan(x2(2, 2).real()));
    CHECK(std::isnan(x2(2, 2).imag()));
}

// TODO: Math ops int

TEST_CASE("su math_ops_double", "[domain]")
{
    using namespace Lielab::domain;
    
    const std::complex<double> j(0.0, 1.0);

    su x1(2);
    x1.set_vector({1.25, 2.5, 3.75});

    const su x1_lm_2 = 2.0*x1;
    CHECK(x1_lm_2(0) == 2.5);
    CHECK(x1_lm_2(1) == 5.0);
    CHECK(x1_lm_2(2) == 7.5);

    const su x1_rm_2 = x1*2.0;
    CHECK(x1_rm_2(0) == 2.5);
    CHECK(x1_rm_2(1) == 5.0);
    CHECK(x1_rm_2(2) == 7.5);

    x1 *= 2.0;
    CHECK(x1(0) == 2.5);
    CHECK(x1(1) == 5.0);
    CHECK(x1(2) == 7.5);

    x1.set_vector({1.25, 2.5, 3.75});
    // TODO: Imaginary ops move elements out of the algebra. Dunno how this should be handled.
    const su x1_lm_2j = (2.0*j)*x1;
    CHECK(x1_lm_2j(0) == 0.0);
    CHECK(x1_lm_2j(1) == 0.0);
    CHECK(x1_lm_2j(2) == 0.0);

    const su x1_rm_2j = x1*(2.0*j);
    CHECK(x1_rm_2j(0) == 0.0);
    CHECK(x1_rm_2j(1) == 0.0);
    CHECK(x1_rm_2j(2) == 0.0);

    x1 *= 2.0*j;
    CHECK(x1(0) == 0.0);
    CHECK(x1(1) == 0.0);
    CHECK(x1(2) == 0.0);

    x1.set_vector({1.25, 2.5, 3.75});

    const su x1_d_2 = x1/2.0;
    CHECK(x1_d_2(0) == 0.625);
    CHECK(x1_d_2(1) == 1.25);
    CHECK(x1_d_2(2) == 1.875);

    x1 /= 2.0;
    CHECK(x1(0) == 0.625);
    CHECK(x1(1) == 1.25);
    CHECK(x1(2) == 1.875);

    x1.set_vector({1.25, 2.5, 3.75});

    const su x1_d_2j = x1/(2.0*j);
    CHECK(x1_d_2j(0) == 0.0);
    CHECK(x1_d_2j(1) == 0.0);
    CHECK(x1_d_2j(2) == 0.0);

    x1 /= 2.0*j;
    CHECK(x1(0) == 0.0);
    CHECK(x1(1) == 0.0);
    CHECK(x1(2) == 0.0);
}

TEST_CASE("su math_ops_cn", "[domain]")
{
    using namespace Lielab::domain;

    su x1(2), x2(2);
    x1.set_vector({1.0, 2.0, 3.0});
    x2.set_vector({1.25, 2.5, 3.75});

    const su x1_add_x2 = x1 + x2;
    CHECK(x1_add_x2(0) == 2.25);
    CHECK(x1_add_x2(1) == 4.5);
    CHECK(x1_add_x2(2) == 6.75);

    x1 += x2;
    CHECK(x1(0) == 2.25);
    CHECK(x1(1) == 4.5);
    CHECK(x1(2) == 6.75);

    x1.set_vector({1.0, 2.0, 3.0});

    const su x1_sub_x2 = x1 - x2;
    CHECK(x1_sub_x2(0) == -0.25);
    CHECK(x1_sub_x2(1) == -0.5);
    CHECK(x1_sub_x2(2) == -0.75);

    x1 -= x2;
    CHECK(x1(0) == -0.25);
    CHECK(x1(1) == -0.5);
    CHECK(x1(2) == -0.75);

    x1.set_vector({1.0, 2.0, 3.0});

    const su x1_unary_sub = (-x1);
    CHECK(x1_unary_sub(0) == -1.0);
    CHECK(x1_unary_sub(1) == -2.0);
    CHECK(x1_unary_sub(2) == -3.0);
}

TEST_CASE("su from_vector", "[domain]")
{
    /*!
    * Tests the from_vector operation.
    */

    using namespace Lielab::domain;

    const su x0 = su::from_vector({});
    const Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0.get_shape() == 1);
    CHECK(x0bar.size() == 0);

    const su x1 = su::from_vector({1.0});
    const Eigen::VectorXd x1bar = x1.get_vector();

    CHECK(x1.get_shape() == 2);
    REQUIRE(x1bar.size() == 3);
    CHECK(x1bar(0) == 1.0);
    CHECK(x1bar(1) == 0.0);
    CHECK(x1bar(2) == 0.0);

    const su x2 = su::from_vector({1.0, 2.0});
    const Eigen::VectorXd x2bar = x2.get_vector();

    CHECK(x2.get_shape() == 2);
    REQUIRE(x2bar.size() == 3);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 0.0);

    const su x3 = su::from_vector({1.0, 2.0, 3.0});
    const Eigen::VectorXd x3bar = x3.get_vector();

    CHECK(x3.get_shape() == 2);
    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
}

// TODO: Project

TEST_CASE("su2", "[domain]")
{
    /*!
    * Tests the su algebra with su(2).
    */

    using namespace Lielab::domain;
    using namespace Lielab::functions;

    const size_t D = su::basis(0, 2).get_dimension();

    // Construct the su2 basis
    std::vector<su> basis;
    for (size_t ii = 0; ii < D; ii++)
    {
        basis.push_back(su::basis(ii, 2));
    }

    // su2 specific identities
    assert_domain(commutator(basis[0], basis[1]),  2*basis[2]);
    assert_domain(commutator(basis[1], basis[2]),  2*basis[0]);
    assert_domain(commutator(basis[2], basis[0]),  2*basis[1]);
    assert_domain(commutator(basis[1], basis[0]), -2*basis[2]);
    assert_domain(commutator(basis[2], basis[1]), -2*basis[0]);
    assert_domain(commutator(basis[0], basis[2]), -2*basis[1]);

    // Hamilton's identities
    // Note that i^2 = j^2 = k^2 = -1^2 isn't checked since this isn't true for the algebra

    // ij = -ji = k
    // assert_domain(b[0]*b[1], -b[1]*b[0]);
    // assert_domain(b[0]*b[1],  b[2]);

    // // jk = -kj = i
    // assert_domain(b[1]*b[2], -b[2]*b[1]);
    // assert_domain(b[1]*b[2],  b[0]);

    // // ki = -ik = j
    // assert_domain(b[2]*b[0], -b[0]*b[2]);
    // assert_domain(b[2]*b[0],  b[1]);
}

TEST_CASE("su3", "[domain]")
{
    /*!
    * Tests the su algebra with su(3).
    */

    using namespace Lielab::domain;

    const size_t D = su::basis(0, 3).get_dimension();

    // Construct the su3 basis
    std::vector<su> b;
    for (size_t ii = 0; ii < D; ii++)
    {
        b.push_back(su::basis(ii, 3));
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
