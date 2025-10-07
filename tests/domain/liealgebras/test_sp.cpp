#include <Lielab.hpp>
#include <iostream>

#include <catch2/catch_all.hpp>

TEST_CASE("sp to_string", "[domain]")
{
    using namespace Lielab::domain;

    const sp x0 = sp(0);
    CHECK(x0.to_string() == "sp(0, R)");
    const sp x1 = sp(1);
    CHECK(x1.to_string() == "sp(0, R)");
    const sp x2 = sp(2);
    CHECK(x2.to_string() == "sp(2, R)");
    const sp x10 = sp(10);
    CHECK(x10.to_string() == "sp(10, R)");
}

TEST_CASE("sp main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const sp xblank = sp();
    CHECK(xblank.get_dimension() == 0);

    const sp x0 = sp(0);
    CHECK(x0.get_dimension() == 0);
    const sp x1 = sp(1);
    CHECK(x1.get_dimension() == 0);
    const sp x3 = sp(3);
    CHECK(x3.get_dimension() == 3);
}

TEST_CASE("sp matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const sp x0 = sp(Eigen::MatrixXd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    CHECK_THROWS(sp(Eigen::MatrixXd::Random(1, 1)));

    const sp x2 = sp(Eigen::MatrixXd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(sp(Eigen::MatrixXd::Random(2, 4)));
    CHECK_THROWS(sp(Eigen::MatrixXd::Random(4, 2)));
}

TEST_CASE("sp basis_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const sp xm10 = sp::basis(-1, 0);
    CHECK(xm10.get_dimension() == 0);
    const Eigen::VectorXd xm10bar = xm10.get_vector();
    CHECK(xm10bar.size() == 0);

    const sp x00 = sp::basis(0, 0);
    CHECK(x00.get_dimension() == 0);
    const Eigen::VectorXd x00bar = x00.get_vector();
    CHECK(x00bar.size() == 0);

    const sp x10 = sp::basis(1, 0);
    CHECK(x10.get_dimension() == 0);
    const Eigen::VectorXd x10bar = x10.get_vector();
    CHECK(x10bar.size() == 0);

    const sp x02 = sp::basis(0, 2);
    CHECK(x02.get_dimension() == 3);
    const Eigen::VectorXd x02bar = x02.get_vector();
    REQUIRE(x02bar.size() == 3);
    CHECK(x02bar(0) == 1.0);
    CHECK(x02bar(1) == 0.0);
    CHECK(x02bar(2) == 0.0);

    const sp x12 = sp::basis(1, 2);
    CHECK(x12.get_dimension() == 3);
    const Eigen::VectorXd x12bar = x12.get_vector();
    REQUIRE(x12bar.size() == 3);
    CHECK(x12bar(0) == 0.0);
    CHECK(x12bar(1) == 1.0);
    CHECK(x12bar(2) == 0.0);

    const sp x22 = sp::basis(2, 2);
    CHECK(x22.get_dimension() == 3);
    const Eigen::VectorXd x22bar = x22.get_vector();
    REQUIRE(x22bar.size() == 3);
    CHECK(x22bar(0) == 0.0);
    CHECK(x22bar(1) == 0.0);
    CHECK(x22bar(2) == 1.0);

    const sp x04 = sp::basis(0, 4);
    CHECK(x04.get_dimension() == 10);
    const Eigen::VectorXd x04bar = x04.get_vector();
    REQUIRE(x04bar.size() == 10);
    CHECK(x04bar(0) == 1.0);
    CHECK(x04bar(1) == 0.0);
    CHECK(x04bar(2) == 0.0);
    CHECK(x04bar(3) == 0.0);
    CHECK(x04bar(4) == 0.0);
    CHECK(x04bar(5) == 0.0);
    CHECK(x04bar(6) == 0.0);
    CHECK(x04bar(7) == 0.0);
    CHECK(x04bar(8) == 0.0);
    CHECK(x04bar(9) == 0.0);
}

TEST_CASE("sp from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const sp x0 = sp::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const sp x1 = sp::from_shape(1);
    CHECK(x1.get_dimension() == 0);
    const Eigen::MatrixXd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 0);
    CHECK(x1hat.cols() == 0);

    const sp x2 = sp::from_shape(2);
    CHECK(x2.get_dimension() == 3);
    const Eigen::MatrixXd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("sp get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    sp zero(0), two(2), four(4), six(6), eight(8);

    // Dimensions
    CHECK(zero.get_dimension() == 0);
    CHECK(two.get_dimension() == 3);
    CHECK(four.get_dimension() == 10);
    CHECK(six.get_dimension() == 21);
    CHECK(eight.get_dimension() == 36);
}

TEST_CASE("sp set/get_vector", "[domain]")
{
    /*!
    * Tests the set/get_vector operation.
    */

    using namespace Lielab::domain;

    sp x0 = sp(0);
    x0.set_vector({});
    Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0bar.size() == 0);

    sp x2 = sp(2);
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

    x2.set_vector({8.0});
    x2bar = x2.get_vector();

    REQUIRE(x2bar.size() == 3);
    CHECK(x2bar(0) == 8.0);
    CHECK(x2bar(1) == 6.0);
    CHECK(x2bar(2) == 7.0);

    sp x4 = sp(4);
    x4.set_vector({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0});
    Eigen::VectorXd x4bar = x4.get_vector();

    REQUIRE(x4bar.size() == 10);
    CHECK(x4bar(0) == 1.0);
    CHECK(x4bar(1) == 2.0);
    CHECK(x4bar(2) == 3.0);
    CHECK(x4bar(3) == 4.0);
    CHECK(x4bar(4) == 5.0);
    CHECK(x4bar(5) == 6.0);
    CHECK(x4bar(6) == 7.0);
    CHECK(x4bar(7) == 8.0);
    CHECK(x4bar(8) == 9.0);
    CHECK(x4bar(9) == 0.0);

    x4.set_vector({10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0});
    x4bar = x4.get_vector();

    REQUIRE(x4bar.size() == 10);
    CHECK(x4bar(0) == 10.0);
    CHECK(x4bar(1) == 11.0);
    CHECK(x4bar(2) == 12.0);
    CHECK(x4bar(3) == 13.0);
    CHECK(x4bar(4) == 14.0);
    CHECK(x4bar(5) == 15.0);
    CHECK(x4bar(6) == 16.0);
    CHECK(x4bar(7) == 17.0);
    CHECK(x4bar(8) == 18.0);
    CHECK(x4bar(9) == 19.0);

    x4.set_vector({20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0});
    x4bar = x4.get_vector();

    REQUIRE(x4bar.size() == 10);
    CHECK(x4bar(0) == 20.0);
    CHECK(x4bar(1) == 21.0);
    CHECK(x4bar(2) == 22.0);
    CHECK(x4bar(3) == 23.0);
    CHECK(x4bar(4) == 24.0);
    CHECK(x4bar(5) == 25.0);
    CHECK(x4bar(6) == 26.0);
    CHECK(x4bar(7) == 27.0);
    CHECK(x4bar(8) == 28.0);
    CHECK(x4bar(9) == 29.0);
}

TEST_CASE("sp get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    sp x0 = sp(0);
    x0.set_vector({});
    Eigen::MatrixXd x0hat = x0.get_matrix();

    REQUIRE(x0hat.rows() == 0);
    REQUIRE(x0hat.cols() == 0);

    sp x2 = sp(2);
    x2.set_vector({1.0, 2.0, 3.0, 4.0});
    Eigen::MatrixXd x2hat = x2.get_matrix();

    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 1.0);
    CHECK(x2hat(0, 1) == 2.0);
    CHECK(x2hat(1, 0) == 3.0);
    CHECK(x2hat(1, 1) == -1.0);

    x2.set_vector({5.0, 6.0, 7.0});
    x2hat = x2.get_matrix();

    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 5.0);
    CHECK(x2hat(0, 1) == 6.0);
    CHECK(x2hat(1, 0) == 7.0);
    CHECK(x2hat(1, 1) == -5.0);

    x2.set_vector({8.0, 9.0});
    x2hat = x2.get_matrix();

    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 8.0);
    CHECK(x2hat(0, 1) == 9.0);
    CHECK(x2hat(1, 0) == 7.0);
    CHECK(x2hat(1, 1) == -8.0);
}

TEST_CASE("sp operator()", "[domain]")
{
    using namespace Lielab::domain;

    sp x0 = sp(0);
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

    sp x2 = sp(2);
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
    CHECK(x2(0, 0) == 1.0);
    CHECK(x2(0, 1) == 2.0);
    CHECK(x2(1, 0) == 3.0);
    CHECK(x2(1, 1) == -1.0);
    CHECK(x2(-1, -1) == -1.0);
    CHECK(x2(-1, -2) == 3.0);
    CHECK(x2(-2, -1) == 2.0);
    CHECK(x2(-2, -2) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x2(0, -3)));
    CHECK(std::isnan(x2(-3, 0)));
    CHECK(std::isnan(x2(-3, -3)));
    CHECK(std::isnan(x2(0, 2)));
    CHECK(std::isnan(x2(2, 0)));
    CHECK(std::isnan(x2(2, 2)));
}

// TODO: Math ops int

TEST_CASE("sp math_ops_double", "[domain]")
{
    using namespace Lielab::domain;

    sp x1(2);
    x1.set_vector({1.25, 2.5, 3.75});

    const sp x1_lm_2 = 2.0*x1;
    CHECK(x1_lm_2(0) == 2.5);
    CHECK(x1_lm_2(1) == 5.0);
    CHECK(x1_lm_2(2) == 7.5);

    const sp x1_rm_2 = x1*2.0;
    CHECK(x1_rm_2(0) == 2.5);
    CHECK(x1_rm_2(1) == 5.0);
    CHECK(x1_rm_2(2) == 7.5);

    x1 *= 2.0;
    CHECK(x1(0) == 2.5);
    CHECK(x1(1) == 5.0);
    CHECK(x1(2) == 7.5);

    x1.set_vector({1.25, 2.5, 3.75});

    const sp x1_d_2 = x1/2.0;
    CHECK(x1_d_2(0) == 0.625);
    CHECK(x1_d_2(1) == 1.25);
    CHECK(x1_d_2(2) == 1.875);

    x1 /= 2.0;
    CHECK(x1(0) == 0.625);
    CHECK(x1(1) == 1.25);
    CHECK(x1(2) == 1.875);
}

TEST_CASE("sp math_ops_sp", "[domain]")
{
    using namespace Lielab::domain;

    sp x1(2), x2(2);
    x1.set_vector({1.0, 2.0, 3.0});
    x2.set_vector({1.25, 2.5, 3.75});

    const sp x1_add_x2 = x1 + x2;
    CHECK(x1_add_x2(0) == 2.25);
    CHECK(x1_add_x2(1) == 4.5);
    CHECK(x1_add_x2(2) == 6.75);

    x1 += x2;
    CHECK(x1(0) == 2.25);
    CHECK(x1(1) == 4.5);
    CHECK(x1(2) == 6.75);

    x1.set_vector({1.0, 2.0, 3.0});

    const sp x1_sub_x2 = x1 - x2;
    CHECK(x1_sub_x2(0) == -0.25);
    CHECK(x1_sub_x2(1) == -0.5);
    CHECK(x1_sub_x2(2) == -0.75);

    x1 -= x2;
    CHECK(x1(0) == -0.25);
    CHECK(x1(1) == -0.5);
    CHECK(x1(2) == -0.75);

    x1.set_vector({1.0, 2.0, 3.0});

    const sp x1_unary_sub = (-x1);
    CHECK(x1_unary_sub(0) == -1.0);
    CHECK(x1_unary_sub(1) == -2.0);
    CHECK(x1_unary_sub(2) == -3.0);
}

TEST_CASE("sp from_vector", "[domain]")
{
    /*!
    * Tests the from_vector operation.
    */

    using namespace Lielab::domain;

    const sp x0 = sp::from_vector({});
    const Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0.get_shape() == 0);
    CHECK(x0bar.size() == 0);

    const sp x1 = sp::from_vector({1.0});
    const Eigen::VectorXd x1bar = x1.get_vector();

    CHECK(x1.get_shape() == 2);
    REQUIRE(x1bar.size() == 3);
    CHECK(x1bar(0) == 1.0);
    CHECK(x1bar(1) == 0.0);
    CHECK(x1bar(2) == 0.0);

    const sp x2 = sp::from_vector({1.0, 2.0});
    const Eigen::VectorXd x2bar = x2.get_vector();

    CHECK(x2.get_shape() == 2);
    REQUIRE(x2bar.size() == 3);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 0.0);

    const sp x3 = sp::from_vector({1.0, 2.0, 3.0});
    const Eigen::VectorXd x3bar = x3.get_vector();

    CHECK(x3.get_shape() == 2);
    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);

    const sp x4 = sp::from_vector({1.0, 2.0, 3.0, 4.0});
    const Eigen::VectorXd x4bar = x4.get_vector();

    CHECK(x4.get_shape() == 4);
    REQUIRE(x4bar.size() == 10);
    CHECK(x4bar(0) == 1.0);
    CHECK(x4bar(1) == 2.0);
    CHECK(x4bar(2) == 3.0);
    CHECK(x4bar(3) == 4.0);
    CHECK(x4bar(4) == 0.0);
    CHECK(x4bar(5) == 0.0);
    CHECK(x4bar(6) == 0.0);
    CHECK(x4bar(7) == 0.0);
    CHECK(x4bar(8) == 0.0);
    CHECK(x4bar(9) == 0.0);
}

TEST_CASE("sp project", "[domain]")
{
    using namespace Lielab::domain;

    const Eigen::MatrixXd rand_2_2 = Eigen::MatrixXd::Random(2, 2);
    const Eigen::MatrixXd proj_2_2 = sp::project(rand_2_2);

    REQUIRE(proj_2_2.rows() == 2);
    REQUIRE(proj_2_2.cols() == 2);
    CHECK(proj_2_2(0, 1) == rand_2_2(0, 1));
    CHECK(proj_2_2(1, 0) == rand_2_2(1, 0));
    CHECK(proj_2_2.trace() == 0.0);

    const Eigen::MatrixXd rand_2_3 = Eigen::MatrixXd::Random(2, 3);
    const Eigen::MatrixXd proj_2_3 = sp::project(rand_2_3);

    REQUIRE(proj_2_3.rows() == 2);
    REQUIRE(proj_2_3.cols() == 2);
    CHECK(proj_2_3(0, 1) == rand_2_3(0, 1));
    CHECK(proj_2_3(1, 0) == rand_2_3(1, 0));
    CHECK(proj_2_3.trace() == 0.0);

    const Eigen::MatrixXd rand_3_2 = Eigen::MatrixXd::Random(3, 2);
    const Eigen::MatrixXd proj_3_2 = sp::project(rand_3_2);

    REQUIRE(proj_3_2.rows() == 2);
    REQUIRE(proj_3_2.cols() == 2);
    CHECK(proj_3_2(0, 1) == rand_3_2(0, 1));
    CHECK(proj_3_2(1, 0) == rand_3_2(1, 0));
    CHECK(proj_3_2.trace() == 0.0);

    const Eigen::MatrixXd rand_3_3 = Eigen::MatrixXd::Random(3, 3);
    const Eigen::MatrixXd proj_3_3 = sp::project(rand_3_3);

    REQUIRE(proj_3_3.rows() == 2);
    REQUIRE(proj_3_3.cols() == 2);
    CHECK(proj_3_3(0, 1) == rand_3_3(0, 1));
    CHECK(proj_3_3(1, 0) == rand_3_3(1, 0));
    CHECK(proj_3_3.trace() == 0.0);
}
