#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

#include <iostream>

#include "../../test_utils.hpp"

TEST_CASE("so to_string", "[domain]")
{
    using namespace Lielab::domain;

    const so x0 = so(0);
    CHECK(x0.to_string() == "so(0)");
    const so x1 = so(1);
    CHECK(x1.to_string() == "so(1)");
    const so x10 = so(10);
    CHECK(x10.to_string() == "so(10)");
}

TEST_CASE("so main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const so xblank = so();
    CHECK(xblank.get_dimension() == 0);

    const so x0 = so(0);
    CHECK(x0.get_dimension() == 0);
    const so x1 = so(1);
    CHECK(x1.get_dimension() == 0);
    const so x3 = so(3);
    CHECK(x3.get_dimension() == 3);
}

TEST_CASE("so matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const so x0 = so(Eigen::MatrixXd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    const so x1 = so(Eigen::MatrixXd::Random(1, 1));
    CHECK(x1.get_shape() == 1);

    const so x2 = so(Eigen::MatrixXd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(so(Eigen::MatrixXd::Random(2, 3)));
    CHECK_THROWS(so(Eigen::MatrixXd::Random(3, 2)));
}

TEST_CASE("so basis_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const so xm10 = so::basis(-1, 0);
    CHECK(xm10.get_dimension() == 0);
    const Eigen::VectorXd xm10bar = xm10.get_vector();
    CHECK(xm10bar.size() == 0);

    const so x00 = so::basis(0, 0);
    CHECK(x00.get_dimension() == 0);
    const Eigen::VectorXd x00bar = x00.get_vector();
    CHECK(x00bar.size() == 0);

    const so x10 = so::basis(1, 0);
    CHECK(x10.get_dimension() == 0);
    const Eigen::VectorXd x10bar = x10.get_vector();
    CHECK(x10bar.size() == 0);

    const so x01 = so::basis(0, 2);
    CHECK(x01.get_dimension() == 1);
    const Eigen::VectorXd x01bar = x01.get_vector();
    REQUIRE(x01bar.size() == 1);
    CHECK(x01bar(0) == 1.0);

    const so x11 = so::basis(1, 2);
    CHECK(x11.get_dimension() == 1);
    const Eigen::VectorXd x11bar = x11.get_vector();
    REQUIRE(x11bar.size() == 1);
    CHECK(x11bar(0) == 0.0);

    const so x21 = so::basis(2, 2);
    CHECK(x21.get_dimension() == 1);
    const Eigen::VectorXd x21bar = x21.get_vector();
    REQUIRE(x21bar.size() == 1);
    CHECK(x21bar(0) == 0.0);

    const so x02 = so::basis(0, 3);
    CHECK(x02.get_dimension() == 3);
    const Eigen::VectorXd x02bar = x02.get_vector();
    REQUIRE(x02bar.size() == 3);
    CHECK(x02bar(0) == 1.0);
    CHECK(x02bar(1) == 0.0);
    CHECK(x02bar(2) == 0.0);
}

TEST_CASE("so from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const so x0 = so::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const so x1 = so::from_shape(1);
    CHECK(x1.get_dimension() == 0);
    const Eigen::MatrixXd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const so x2 = so::from_shape(2);
    CHECK(x2.get_dimension() == 1);
    const Eigen::MatrixXd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("so get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    so zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(zero.get_dimension() == 0);
    CHECK(one.get_dimension() == 0);
    CHECK(two.get_dimension() == 1);
    CHECK(three.get_dimension() == 3);
    CHECK(four.get_dimension() == 6);
    CHECK(five.get_dimension() == 10);
    CHECK(six.get_dimension() == 15);
    CHECK(seven.get_dimension() == 21);
    CHECK(eight.get_dimension() == 28);
}

TEST_CASE("so set/get_vector", "[domain]")
{
    /*!
    * Tests the set/get_vector operation.
    */

    using namespace Lielab::domain;

    so x0 = so(0);
    x0.set_vector({});
    Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0bar.size() == 0);

    so x3 = so(3);
    x3.set_vector({1.0, 2.0, 3.0, 4.0});
    Eigen::VectorXd x3bar = x3.get_vector();

    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);

    x3.set_vector({5.0, 6.0, 7.0});
    x3bar = x3.get_vector();

    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 5.0);
    CHECK(x3bar(1) == 6.0);
    CHECK(x3bar(2) == 7.0);

    x3.set_vector({8.0});
    x3bar = x3.get_vector();

    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 8.0);
    CHECK(x3bar(1) == 6.0);
    CHECK(x3bar(2) == 7.0);

    so x4 = so(4);
    x4.set_vector({1.0, 2.0, 3.0, 4.0, 5.0});
    Eigen::VectorXd x4bar = x4.get_vector();

    REQUIRE(x4bar.size() == 6);
    CHECK(x4bar(0) == 1.0);
    CHECK(x4bar(1) == 2.0);
    CHECK(x4bar(2) == 3.0);
    CHECK(x4bar(3) == 4.0);
    CHECK(x4bar(4) == 5.0);
    CHECK(x4bar(5) == 0.0);

    x4.set_vector({7.0, 8.0, 9.0, 10.0, 11.0, 12.0});
    x4bar = x4.get_vector();

    REQUIRE(x4bar.size() == 6);
    CHECK(x4bar(0) == 7.0);
    CHECK(x4bar(1) == 8.0);
    CHECK(x4bar(2) == 9.0);
    CHECK(x4bar(3) == 10.0);
    CHECK(x4bar(4) == 11.0);
    CHECK(x4bar(5) == 12.0);

    x4.set_vector({13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0});
    x4bar = x4.get_vector();

    REQUIRE(x4bar.size() == 6);
    CHECK(x4bar(0) == 13.0);
    CHECK(x4bar(1) == 14.0);
    CHECK(x4bar(2) == 15.0);
    CHECK(x4bar(3) == 16.0);
    CHECK(x4bar(4) == 17.0);
    CHECK(x4bar(5) == 18.0);
}

TEST_CASE("so get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    so x0 = so(0);
    x0.set_vector({});
    Eigen::MatrixXd x0hat = x0.get_matrix();

    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    so x1 = so(1);
    x1.set_vector({});
    Eigen::MatrixXd x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == 0.0);

    so x2 = so(2);
    x2.set_vector({1.0, 2.0});
    Eigen::MatrixXd x2hat = x2.get_matrix();

    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 0.0);
    CHECK(x2hat(0, 1) == -1.0);
    CHECK(x2hat(1, 0) == 1.0);
    CHECK(x2hat(1, 1) == 0.0);

    x2.set_vector({3.0});
    x2hat = x2.get_matrix();

    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 0.0);
    CHECK(x2hat(0, 1) == -3.0);
    CHECK(x2hat(1, 0) == 3.0);
    CHECK(x2hat(1, 1) == 0.0);

    x2.set_vector({});
    x2hat = x2.get_matrix();

    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 0.0);
    CHECK(x2hat(0, 1) == -3.0);
    CHECK(x2hat(1, 0) == 3.0);
    CHECK(x2hat(1, 1) == 0.0);

    so x3 = so(3);
    x3.set_vector({1.0});
    Eigen::MatrixXd x3hat = x3.get_matrix();
    REQUIRE(x3hat.rows() == 3);
    REQUIRE(x3hat.cols() == 3);
    CHECK(x3hat(0, 0) == 0.0);
    CHECK(x3hat(0, 1) == 0.0);
    CHECK(x3hat(0, 2) == 0.0);
    CHECK(x3hat(1, 0) == 0.0);
    CHECK(x3hat(1, 1) == 0.0);
    CHECK(x3hat(1, 2) == -1.0);
    CHECK(x3hat(2, 0) == 0.0);
    CHECK(x3hat(2, 1) == 1.0);
    CHECK(x3hat(2, 2) == 0.0);

    x3.set_vector({2.0, 3.0});
    x3hat = x3.get_matrix();
    REQUIRE(x3hat.rows() == 3);
    REQUIRE(x3hat.cols() == 3);
    CHECK(x3hat(0, 0) == 0.0);
    CHECK(x3hat(0, 1) == 0.0);
    CHECK(x3hat(0, 2) == 3.0);
    CHECK(x3hat(1, 0) == 0.0);
    CHECK(x3hat(1, 1) == 0.0);
    CHECK(x3hat(1, 2) == -2.0);
    CHECK(x3hat(2, 0) == -3.0);
    CHECK(x3hat(2, 1) == 2.0);
    CHECK(x3hat(2, 2) == 0.0);

    x3.set_vector({4.0, 5.0, 6.0});
    x3hat = x3.get_matrix();
    REQUIRE(x3hat.rows() == 3);
    REQUIRE(x3hat.cols() == 3);
    CHECK(x3hat(0, 0) == 0.0);
    CHECK(x3hat(0, 1) == -6.0);
    CHECK(x3hat(0, 2) == 5.0);
    CHECK(x3hat(1, 0) == 6.0);
    CHECK(x3hat(1, 1) == 0.0);
    CHECK(x3hat(1, 2) == -4.0);
    CHECK(x3hat(2, 0) == -5.0);
    CHECK(x3hat(2, 1) == 4.0);
    CHECK(x3hat(2, 2) == 0.0);
}

TEST_CASE("so operator()", "[domain]")
{
    using namespace Lielab::domain;

    so x0 = so(0);
    x0.set_vector({});

    // Out of bounds
    CHECK(std::isnan(x0(-1)));
    CHECK(std::isnan(x0(0)));
    CHECK(std::isnan(x0(1)));

    // Out of bounds
    CHECK(std::isnan(x0(0, -1)));
    CHECK(std::isnan(x0(-1, 0)));
    CHECK(std::isnan(x0(-1, -1)));
    CHECK(std::isnan(x0(0, 0)));
    CHECK(std::isnan(x0(0, 1)));
    CHECK(std::isnan(x0(1, 0)));
    CHECK(std::isnan(x0(1, 1)));

    so x1 = so(1);
    x1.set_vector({});

    // Out of bounds
    CHECK(std::isnan(x1(-1)));
    CHECK(std::isnan(x1(0)));
    CHECK(std::isnan(x1(1)));

    // In bounds
    CHECK(x1(0, 0) == 0.0);
    CHECK(x1(-1, -1) == 0.0);

    // Out of bounds
    CHECK(std::isnan(x1(0, -2)));
    CHECK(std::isnan(x1(-2, 0)));
    CHECK(std::isnan(x1(-2, -2)));
    CHECK(std::isnan(x1(0, 1)));
    CHECK(std::isnan(x1(1, 0)));
    CHECK(std::isnan(x1(1, 1)));

    so x3 = so(3);
    x3.set_vector({1.0, 2.0, 3.0});

    // In bounds
    CHECK(x3(0) == 1.0);
    CHECK(x3(1) == 2.0);
    CHECK(x3(2) == 3.0);
    CHECK(x3(-1) == 3.0);
    CHECK(x3(-2) == 2.0);
    CHECK(x3(-3) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x3(-4)));
    CHECK(std::isnan(x3(3)));

    // In bounds
    CHECK(x3(0, 0) == 0.0);
    CHECK(x3(0, 1) == -3.0);
    CHECK(x3(0, 2) == 2.0);
    CHECK(x3(1, 0) == 3.0);
    CHECK(x3(1, 1) == 0.0);
    CHECK(x3(1, 2) == -1.0);
    CHECK(x3(2, 0) == -2.0);
    CHECK(x3(2, 1) == 1.0);
    CHECK(x3(2, 2) == 0.0);
    CHECK(x3(-1, -1) == 0.0);
    CHECK(x3(-1, -2) == 1.0);
    CHECK(x3(-1, -3) == -2.0);
    CHECK(x3(-2, -1) == -1.0);
    CHECK(x3(-2, -2) == 0.0);
    CHECK(x3(-2, -3) == 3.0);
    CHECK(x3(-3, -1) == 2.0);
    CHECK(x3(-3, -2) == -3.0);
    CHECK(x3(-3, -3) == 0.0);

    // Out of bounds
    CHECK(std::isnan(x3(0, -4)));
    CHECK(std::isnan(x3(-4, 0)));
    CHECK(std::isnan(x3(-4, -4)));
    CHECK(std::isnan(x3(0, 3)));
    CHECK(std::isnan(x3(3, 0)));
    CHECK(std::isnan(x3(3, 3)));
}

// TODO: Math ops int

TEST_CASE("so math_ops_double", "[domain]")
{
    using namespace Lielab::domain;

    so x1(3);
    x1.set_vector({1.25, 2.5, 3.75});

    const so x1_lm_2 = 2.0*x1;
    CHECK(x1_lm_2(0) == 2.5);
    CHECK(x1_lm_2(1) == 5.0);
    CHECK(x1_lm_2(2) == 7.5);

    const so x1_rm_2 = x1*2.0;
    CHECK(x1_rm_2(0) == 2.5);
    CHECK(x1_rm_2(1) == 5.0);
    CHECK(x1_rm_2(2) == 7.5);

    x1 *= 2.0;
    CHECK(x1(0) == 2.5);
    CHECK(x1(1) == 5.0);
    CHECK(x1(2) == 7.5);

    x1.set_vector({1.25, 2.5, 3.75});

    const so x1_d_2 = x1/2.0;
    CHECK(x1_d_2(0) == 0.625);
    CHECK(x1_d_2(1) == 1.25);
    CHECK(x1_d_2(2) == 1.875);

    x1 /= 2.0;
    CHECK(x1(0) == 0.625);
    CHECK(x1(1) == 1.25);
    CHECK(x1(2) == 1.875);
}

TEST_CASE("so math_ops_so", "[domain]")
{
    using namespace Lielab::domain;

    so x1(3), x2(3);
    x1.set_vector({1.0, 2.0, 3.0});
    x2.set_vector({1.25, 2.5, 3.75});

    const so x1_add_x2 = x1 + x2;
    CHECK(x1_add_x2(0) == 2.25);
    CHECK(x1_add_x2(1) == 4.5);
    CHECK(x1_add_x2(2) == 6.75);

    x1 += x2;
    CHECK(x1(0) == 2.25);
    CHECK(x1(1) == 4.5);
    CHECK(x1(2) == 6.75);

    x1.set_vector({1.0, 2.0, 3.0});

    const so x1_sub_x2 = x1 - x2;
    CHECK(x1_sub_x2(0) == -0.25);
    CHECK(x1_sub_x2(1) == -0.5);
    CHECK(x1_sub_x2(2) == -0.75);

    x1 -= x2;
    CHECK(x1(0) == -0.25);
    CHECK(x1(1) == -0.5);
    CHECK(x1(2) == -0.75);

    x1.set_vector({1.0, 2.0, 3.0});

    const so x1_unary_sub = (-x1);
    CHECK(x1_unary_sub(0) == -1.0);
    CHECK(x1_unary_sub(1) == -2.0);
    CHECK(x1_unary_sub(2) == -3.0);
}

TEST_CASE("so from_vector", "[domain]")
{
    /*!
    * Tests the from_vector operation.
    */

    using namespace Lielab::domain;

    const so x0 = so::from_vector({});
    const Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0.get_shape() == 1);
    CHECK(x0bar.size() == 0);

    const so x1 = so::from_vector({1.0});
    const Eigen::VectorXd x1bar = x1.get_vector();

    CHECK(x1.get_shape() == 2);
    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == 1.0);

    const so x2 = so::from_vector({1.0, 2.0});
    const Eigen::VectorXd x2bar = x2.get_vector();

    CHECK(x2.get_shape() == 3);
    REQUIRE(x2bar.size() == 3);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 0.0);

    const so x3 = so::from_vector({1.0, 2.0, 3.0});
    const Eigen::VectorXd x3bar = x3.get_vector();

    CHECK(x3.get_shape() == 3);
    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
}

TEST_CASE("so project", "[domain]")
{
    using namespace Lielab::domain;

    const Eigen::MatrixXd rand_2_2 = Eigen::MatrixXd::Random(2, 2);
    const Eigen::MatrixXd proj_2_2 = so::project(rand_2_2);

    REQUIRE(proj_2_2.rows() == 2);
    REQUIRE(proj_2_2.cols() == 2);
    CHECK(proj_2_2(0, 0) == 0.0);
    CHECK(proj_2_2(0, 1) == -proj_2_2(1, 0));
    CHECK(proj_2_2(1, 0) == -proj_2_2(0, 1));
    CHECK(proj_2_2(1, 1) == 0.0);

    const Eigen::MatrixXd rand_3_3 = Eigen::MatrixXd::Random(3, 3);
    const Eigen::MatrixXd proj_3_3 = so::project(rand_3_3);

    REQUIRE(proj_3_3.rows() == 3);
    REQUIRE(proj_3_3.cols() == 3);
    CHECK(proj_3_3(0, 0) == 0.0);
    CHECK(proj_3_3(0, 1) == -proj_3_3(1, 0));
    CHECK(proj_3_3(0, 2) == -proj_3_3(2, 0));
    CHECK(proj_3_3(1, 0) == -proj_3_3(0, 1));
    CHECK(proj_3_3(1, 1) == 0.0);
    CHECK(proj_3_3(1, 2) == -proj_3_3(2, 1));
    CHECK(proj_3_3(2, 0) == -proj_3_3(0, 2));
    CHECK(proj_3_3(2, 1) == -proj_3_3(1, 2));
    CHECK(proj_3_3(2, 2) == 0.0);

    const Eigen::MatrixXd rand_2_3 = Eigen::MatrixXd::Random(2, 3);
    const Eigen::MatrixXd proj_2_3 = so::project(rand_2_3);

    REQUIRE(proj_2_3.rows() == 2);
    REQUIRE(proj_2_3.cols() == 2);
    CHECK(proj_2_3(0, 0) == 0.0);
    CHECK(proj_2_3(0, 1) == -proj_2_3(1, 0));
    CHECK(proj_2_3(1, 0) == -proj_2_3(0, 1));
    CHECK(proj_2_3(1, 1) == 0.0);

    const Eigen::MatrixXd rand_3_2 = Eigen::MatrixXd::Random(3, 2);
    const Eigen::MatrixXd proj_3_2 = so::project(rand_3_2);

    REQUIRE(proj_3_2.rows() == 2);
    REQUIRE(proj_3_2.cols() == 2);
    CHECK(proj_3_2(0, 0) == 0.0);
    CHECK(proj_3_2(0, 1) == -proj_3_2(1, 0));
    CHECK(proj_3_2(1, 0) == -proj_3_2(0, 1));
    CHECK(proj_3_2(1, 1) == 0.0);
}

TEST_CASE("so2", "[domain]")
{
    /*!
    Tests the so algebra with so(2).
    
    Note that some of these test cases are trivial since so(2) is 1-dimensional.
    */
    
    using namespace Lielab::domain;
    using namespace Lielab::functions;

    so x = so::basis(0,2);
    so zero = x*0;

    assert_domain(commutator(x, x), zero);
}

TEST_CASE("so3", "[domain]")
{
    /*!
    * Tests the so algebra with so(3).
    */

    using namespace Lielab::domain;
    using namespace Lielab::functions;

    so x = so::basis(0,3);
    so y = so::basis(1,3);
    so z = so::basis(2,3);
    so zero = x*0;

    assert_domain(commutator(x, y), z);
    assert_domain(commutator(y, z), x);
    assert_domain(commutator(z, x), y);
    assert_domain(commutator(y, x), -z);
    assert_domain(commutator(z, y), -x);
    assert_domain(commutator(x, z), -y);
}
