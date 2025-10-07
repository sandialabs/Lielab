#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

#include <iostream>

TEST_CASE("glr to_string", "[domain]")
{
    using namespace Lielab::domain;

    const glr x0 = glr(0);
    CHECK(x0.to_string() == "gl(0, R)");
    const glr x1 = glr(1);
    CHECK(x1.to_string() == "gl(1, R)");
    const glr x10 = glr(10);
    CHECK(x10.to_string() == "gl(10, R)");
}

TEST_CASE("glr main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const glr xblank = glr();
    CHECK(xblank.get_dimension() == 0);

    const glr x0 = glr(0);
    CHECK(x0.get_dimension() == 0);
    const glr x1 = glr(1);
    CHECK(x1.get_dimension() == 1);
    const glr x10 = glr(10);
    CHECK(x10.get_dimension() == 100);
}

TEST_CASE("glr matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const glr x0 = glr(Eigen::MatrixXd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    const glr x1 = glr(Eigen::MatrixXd::Random(1, 1));
    CHECK(x1.get_shape() == 1);

    const glr x2 = glr(Eigen::MatrixXd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(glr(Eigen::MatrixXd::Random(2, 3)));
    CHECK_THROWS(glr(Eigen::MatrixXd::Random(3, 2)));
}

TEST_CASE("glr basis_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const glr xm10 = glr::basis(-1, 0);
    CHECK(xm10.get_dimension() == 0);
    const Eigen::VectorXd xm10bar = xm10.get_vector();
    CHECK(xm10bar.size() == 0);

    const glr x00 = glr::basis(0, 0);
    CHECK(x00.get_dimension() == 0);
    const Eigen::VectorXd x00bar = x00.get_vector();
    CHECK(x00bar.size() == 0);

    const glr x10 = glr::basis(1, 0);
    CHECK(x10.get_dimension() == 0);
    const Eigen::VectorXd x10bar = x10.get_vector();
    CHECK(x10bar.size() == 0);

    const glr x01 = glr::basis(0, 1);
    CHECK(x01.get_dimension() == 1);
    const Eigen::VectorXd x01bar = x01.get_vector();
    REQUIRE(x01bar.size() == 1);
    CHECK(x01bar(0) == 1.0);

    const glr x11 = glr::basis(1, 1);
    CHECK(x11.get_dimension() == 1);
    const Eigen::VectorXd x11bar = x11.get_vector();
    REQUIRE(x11bar.size() == 1);
    CHECK(x11bar(0) == 0.0);

    const glr x02 = glr::basis(0, 2);
    CHECK(x02.get_dimension() == 4);
    const Eigen::VectorXd x02bar = x02.get_vector();
    REQUIRE(x02bar.size() == 4);
    CHECK(x02bar(0) == 1.0);
    CHECK(x02bar(1) == 0.0);
    CHECK(x02bar(2) == 0.0);
    CHECK(x02bar(3) == 0.0);
}

TEST_CASE("glr from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const glr x0 = glr::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const glr x1 = glr::from_shape(1);
    CHECK(x1.get_dimension() == 1);
    const Eigen::MatrixXd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const glr x2 = glr::from_shape(2);
    CHECK(x2.get_dimension() == 4);
    const Eigen::MatrixXd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("glr get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    glr zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(zero.get_dimension() == 0);
    CHECK(one.get_dimension() == 1);
    CHECK(two.get_dimension() == 4);
    CHECK(three.get_dimension() == 9);
    CHECK(four.get_dimension() == 16);
    CHECK(five.get_dimension() == 25);
    CHECK(six.get_dimension() == 36);
    CHECK(seven.get_dimension() == 49);
    CHECK(eight.get_dimension() == 64);
}

TEST_CASE("glr set/get_vector", "[domain]")
{
    /*!
    * Tests the set/get_vector operation.
    */

    using namespace Lielab::domain;

    glr x0 = glr(0);
    x0.set_vector({});
    Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0bar.size() == 0);

    glr x1 = glr(1);
    x1.set_vector({1.0, 2.0});
    Eigen::VectorXd x1bar = x1.get_vector();

    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == 1.0);

    x1.set_vector({3.0});
    x1bar = x1.get_vector();

    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == 3.0);

    x1.set_vector({});
    x1bar = x1.get_vector();

    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == 3.0);

    glr x2 = glr(2);
    x2.set_vector({1.0, 2.0, 3.0});
    Eigen::VectorXd x2bar = x2.get_vector();

    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 3.0);
    CHECK(x2bar(3) == 0.0);

    x2.set_vector({4.0, 5.0, 6.0, 7.0});
    x2bar = x2.get_vector();

    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 4.0);
    CHECK(x2bar(1) == 5.0);
    CHECK(x2bar(2) == 6.0);
    CHECK(x2bar(3) == 7.0);

    x2.set_vector({8.0, 9.0, 10.0, 11.0, 12.0});
    x2bar = x2.get_vector();

    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 8.0);
    CHECK(x2bar(1) == 9.0);
    CHECK(x2bar(2) == 10.0);
    CHECK(x2bar(3) == 11.0);
}

TEST_CASE("glr get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    glr x0 = glr(0);
    x0.set_vector({});
    Eigen::MatrixXd x0hat = x0.get_matrix();

    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    glr x1 = glr(1);
    x1.set_vector({});
    Eigen::MatrixXd x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == 0.0);

    x1.set_vector({1.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == 1.0);

    x1.set_vector({2.0, 3.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == 2.0);

    glr x2 = glr(2);
    x2.set_vector({1.0, 2.0, 3.0});
    Eigen::MatrixXd x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 1.0);
    CHECK(x2hat(0, 1) == 2.0);
    CHECK(x2hat(1, 0) == 3.0);
    CHECK(x2hat(1, 1) == 0.0);

    x2.set_vector({4.0, 5.0, 6.0, 7.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 4.0);
    CHECK(x2hat(0, 1) == 5.0);
    CHECK(x2hat(1, 0) == 6.0);
    CHECK(x2hat(1, 1) == 7.0);

    x2.set_vector({8.0, 9.0, 10.0, 11.0, 12.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 8.0);
    CHECK(x2hat(0, 1) == 9.0);
    CHECK(x2hat(1, 0) == 10.0);
    CHECK(x2hat(1, 1) == 11.0);
}

TEST_CASE("glr operator()", "[domain]")
{
    using namespace Lielab::domain;

    glr x0 = glr::from_shape(0);
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

    glr x1 = glr(1);
    x1.set_vector({1.0});

    // In bounds
    CHECK(x1(0) == 1.0);
    CHECK(x1(-1) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x1(-2)));
    CHECK(std::isnan(x1(1)));

    // In bounds
    CHECK(x1(0, 0) == 1.0);
    CHECK(x1(-1, -1) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x1(0, -2)));
    CHECK(std::isnan(x1(-2, 0)));
    CHECK(std::isnan(x1(-2, -2)));
    CHECK(std::isnan(x1(0, 2)));
    CHECK(std::isnan(x1(2, 0)));
    CHECK(std::isnan(x1(2, 2)));

    glr x2 = glr(2);
    x2.set_vector({1.0, 2.0, 3.0, 4.0});

    // In bounds
    CHECK(x2(0) == 1.0);
    CHECK(x2(1) == 2.0);
    CHECK(x2(2) == 3.0);
    CHECK(x2(3) == 4.0);
    CHECK(x2(-1) == 4.0);
    CHECK(x2(-2) == 3.0);
    CHECK(x2(-3) == 2.0);
    CHECK(x2(-4) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x2(-5)));
    CHECK(std::isnan(x2(4)));

    // In bounds
    CHECK(x2(0, 0) == 1.0);
    CHECK(x2(0, 1) == 2.0);
    CHECK(x2(1, 0) == 3.0);
    CHECK(x2(1, 1) == 4.0);
    CHECK(x2(-1, -1) == 4.0);
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

// TODO: math ops int

TEST_CASE("glr math_ops_double", "[domain]")
{
    using namespace Lielab::domain;

    glr x1(2);
    x1.set_vector({1.25, 2.5, 3.75, 1.0});

    const glr x1_lm_2 = 2.0*x1;
    CHECK(x1_lm_2(0, 0) == 2.5);
    CHECK(x1_lm_2(0, 1) == 5.0);
    CHECK(x1_lm_2(1, 0) == 7.5);
    CHECK(x1_lm_2(1, 1) == 2.0);

    const glr x1_rm_2 = x1*2.0;
    CHECK(x1_rm_2(0, 0) == 2.5);
    CHECK(x1_rm_2(0, 1) == 5.0);
    CHECK(x1_rm_2(1, 0) == 7.5);
    CHECK(x1_rm_2(1, 1) == 2.0);

    x1 *= 2.0;
    CHECK(x1(0, 0) == 2.5);
    CHECK(x1(0, 1) == 5.0);
    CHECK(x1(1, 0) == 7.5);
    CHECK(x1(1, 1) == 2.0);

    x1.set_vector({1.25, 2.5, 3.75, 1.0});

    const glr x1_d_2 = x1/2.0;
    CHECK(x1_d_2(0, 0) == 0.625);
    CHECK(x1_d_2(0, 1) == 1.25);
    CHECK(x1_d_2(1, 0) == 1.875);
    CHECK(x1_d_2(1, 1) == 0.5);

    x1 /= 2.0;
    CHECK(x1(0, 0) == 0.625);
    CHECK(x1(0, 1) == 1.25);
    CHECK(x1(1, 0) == 1.875);
    CHECK(x1(1, 1) == 0.5);
}

TEST_CASE("glr math_ops_cn", "[domain]")
{
    using namespace Lielab::domain;

    glr x1(2), x2(2);
    x1.set_vector({1.0, 2.0, 3.0, 4.0});
    x2.set_vector({1.25, 2.5, 3.75, 1.0});

    const glr x1_add_x2 = x1 + x2;
    CHECK(x1_add_x2(0, 0) == 2.25);
    CHECK(x1_add_x2(0, 1) == 4.5);
    CHECK(x1_add_x2(1, 0) == 6.75);
    CHECK(x1_add_x2(1, 1) == 5.0);

    x1 += x2;
    CHECK(x1(0, 0) == 2.25);
    CHECK(x1(0, 1) == 4.5);
    CHECK(x1(1, 0) == 6.75);
    CHECK(x1(1, 1) == 5.0);

    x1.set_vector({1.0, 2.0, 3.0, 4.0});

    const glr x1_sub_x2 = x1 - x2;
    CHECK(x1_sub_x2(0, 0) == -0.25);
    CHECK(x1_sub_x2(0, 1) == -0.5);
    CHECK(x1_sub_x2(1, 0) == -0.75);
    CHECK(x1_sub_x2(1, 1) == 3.0);

    x1 -= x2;
    CHECK(x1(0, 0) == -0.25);
    CHECK(x1(0, 1) == -0.5);
    CHECK(x1(1, 0) == -0.75);
    CHECK(x1(1, 1) == 3.0);

    x1.set_vector({1.0, 2.0, 3.0, 4.0});

    const glr x1_unary_sub = (-x1);
    CHECK(x1_unary_sub(0, 0) == -1.0);
    CHECK(x1_unary_sub(0, 1) == -2.0);
    CHECK(x1_unary_sub(1, 0) == -3.0);
    CHECK(x1_unary_sub(1, 1) == -4.0);
}

TEST_CASE("glr get/from_vector", "[domain]")
{
    /*!
    * Tests the get/from_vector operation.
    */

    using namespace Lielab::domain;

    const glr x0 = glr::from_vector({});
    const Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0.get_shape() == 0);
    CHECK(x0bar.size() == 0);

    const glr x1 = glr::from_vector({1.0});
    const Eigen::VectorXd x1bar = x1.get_vector();

    CHECK(x1.get_shape() == 1);
    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == 1.0);

    const glr x2 = glr::from_vector({1.0, 2.0});
    const Eigen::VectorXd x2bar = x2.get_vector();

    CHECK(x2.get_shape() == 2);
    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 0.0);
    CHECK(x2bar(3) == 0.0);

    const glr x3 = glr::from_vector({1.0, 2.0, 3.0});
    const Eigen::VectorXd x3bar = x3.get_vector();

    CHECK(x3.get_shape() == 2);
    REQUIRE(x3bar.size() == 4);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
    CHECK(x3bar(3) == 0.0);

    const glr x5 = glr::from_vector({1.0, 2.0, 3.0, 4.0, 5.0});
    const Eigen::VectorXd x5bar = x5.get_vector();

    CHECK(x5.get_shape() == 3);
    REQUIRE(x5bar.size() == 9);
    CHECK(x5bar(0) == 1.0);
    CHECK(x5bar(1) == 2.0);
    CHECK(x5bar(2) == 3.0);
    CHECK(x5bar(3) == 4.0);
    CHECK(x5bar(4) == 5.0);
    CHECK(x5bar(5) == 0.0);
    CHECK(x5bar(6) == 0.0);
    CHECK(x5bar(7) == 0.0);
    CHECK(x5bar(8) == 0.0);
}

TEST_CASE("glr project", "[domain]")
{
    using namespace Lielab::domain;

    const Eigen::MatrixXd rand_2_2 = Eigen::MatrixXd::Random(2, 2);
    const Eigen::MatrixXd proj_2_2 = glr::project(rand_2_2);

    REQUIRE(proj_2_2.rows() == 2);
    REQUIRE(proj_2_2.cols() == 2);
    CHECK(proj_2_2(0, 0) == rand_2_2(0, 0));
    CHECK(proj_2_2(0, 1) == rand_2_2(0, 1));
    CHECK(proj_2_2(1, 0) == rand_2_2(1, 0));
    CHECK(proj_2_2(1, 1) == rand_2_2(1, 1));

    const Eigen::MatrixXd rand_3_3 = Eigen::MatrixXd::Random(3, 3);
    const Eigen::MatrixXd proj_3_3 = glr::project(rand_3_3);

    REQUIRE(proj_3_3.rows() == 3);
    REQUIRE(proj_3_3.cols() == 3);
    CHECK(proj_3_3(0, 0) == rand_3_3(0, 0));
    CHECK(proj_3_3(0, 1) == rand_3_3(0, 1));
    CHECK(proj_3_3(0, 2) == rand_3_3(0, 2));
    CHECK(proj_3_3(1, 0) == rand_3_3(1, 0));
    CHECK(proj_3_3(1, 1) == rand_3_3(1, 1));
    CHECK(proj_3_3(1, 2) == rand_3_3(1, 2));
    CHECK(proj_3_3(2, 0) == rand_3_3(2, 0));
    CHECK(proj_3_3(2, 1) == rand_3_3(2, 1));
    CHECK(proj_3_3(2, 2) == rand_3_3(2, 2));

    const Eigen::MatrixXd rand_2_3 = Eigen::MatrixXd::Random(2, 3);
    const Eigen::MatrixXd proj_2_3 = glr::project(rand_2_3);

    REQUIRE(proj_2_3.rows() == 2);
    REQUIRE(proj_2_3.cols() == 2);
    CHECK(proj_2_3(0, 0) == rand_2_3(0, 0));
    CHECK(proj_2_3(0, 1) == rand_2_3(0, 1));
    CHECK(proj_2_3(1, 0) == rand_2_3(1, 0));
    CHECK(proj_2_3(1, 1) == rand_2_3(1, 1));

    const Eigen::MatrixXd rand_3_2 = Eigen::MatrixXd::Random(3, 2);
    const Eigen::MatrixXd proj_3_2 = glr::project(rand_3_2);

    REQUIRE(proj_3_2.rows() == 2);
    REQUIRE(proj_3_2.cols() == 2);
    CHECK(proj_3_2(0, 0) == rand_3_2(0, 0));
    CHECK(proj_3_2(0, 1) == rand_3_2(0, 1));
    CHECK(proj_3_2(1, 0) == rand_3_2(1, 0));
    CHECK(proj_3_2(1, 1) == rand_3_2(1, 1));
}
