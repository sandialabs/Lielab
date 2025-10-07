#include <Lielab.hpp>
#include <iostream>

#include <catch2/catch_all.hpp>

TEST_CASE("RN to_string", "[domain]")
{
    using namespace Lielab::domain;

    const RN xzero = RN::from_shape(0);
    CHECK(xzero.to_string() == "R^nan");
    const RN x0 = RN(0);
    CHECK(x0.to_string() == "R^0");
    const RN x1 = RN(1);
    CHECK(x1.to_string() == "R^1");
    const RN x10 = RN(10);
    CHECK(x10.to_string() == "R^10");
}

TEST_CASE("RN main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const RN xblank = RN();
    CHECK(xblank.get_dimension() == 0);

    const RN x0 = RN(0);
    CHECK(x0.get_dimension() == 0);
    const RN x1 = RN(1);
    CHECK(x1.get_dimension() == 1);
    const RN x10 = RN(10);
    CHECK(x10.get_dimension() == 10);
}

TEST_CASE("RN matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const RN x0 = RN(Eigen::MatrixXd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    const RN x1 = RN(Eigen::MatrixXd::Random(1, 1));
    CHECK(x1.get_shape() == 1);

    const RN x2 = RN(Eigen::MatrixXd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(RN(Eigen::MatrixXd::Random(2, 3)));
    CHECK_THROWS(RN(Eigen::MatrixXd::Random(3, 2)));
}

TEST_CASE("RN from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const RN x0 = RN::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const RN x1 = RN::from_shape(1);
    CHECK(x1.get_dimension() == 0);
    const Eigen::MatrixXd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const RN x2 = RN::from_shape(2);
    CHECK(x2.get_dimension() == 1);
    const Eigen::MatrixXd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("RN get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    RN veryzero = RN::from_shape(0);
    RN zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    CHECK(veryzero.get_dimension() == 0);
    CHECK(zero.get_dimension() == 0);
    CHECK(one.get_dimension() == 1);
    CHECK(two.get_dimension() == 2);
    CHECK(three.get_dimension() == 3);
    CHECK(four.get_dimension() == 4);
    CHECK(five.get_dimension() == 5);
    CHECK(six.get_dimension() == 6);
    CHECK(seven.get_dimension() == 7);
    CHECK(eight.get_dimension() == 8);
}

TEST_CASE("RN serialize/unserialize", "[domain]")
{
    /*!
    * Tests the serialize/unserialize operation.
    */

    using namespace Lielab::domain;

    RN xzero = RN::from_shape(0);
    xzero.unserialize({});
    Eigen::VectorXd xzerobar = xzero.serialize();

    CHECK(xzerobar.size() == 0);

    RN x0 = RN(0);
    x0.unserialize({});
    Eigen::VectorXd x0bar = x0.serialize();

    CHECK(x0bar.size() == 0);

    RN x2 = RN(2);
    x2.unserialize({1.0, 2.0, 3.0});
    Eigen::VectorXd x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 2);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);

    x2.unserialize({4.0, 5.0});
    x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 2);
    CHECK(x2bar(0) == 4.0);
    CHECK(x2bar(1) == 5.0);

    x2.unserialize({6.0});
    x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 2);
    CHECK(x2bar(0) == 6.0);
    CHECK(x2bar(1) == 5.0);

    RN x3 = RN(3);
    x3.unserialize({1.0, 2.0});
    Eigen::VectorXd x3bar = x3.serialize();

    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 0.0);

    x3.unserialize({3.0, 4.0, 5.0});
    x3bar = x3.serialize();

    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 3.0);
    CHECK(x3bar(1) == 4.0);
    CHECK(x3bar(2) == 5.0);

    x3.unserialize({6.0, 7.0, 8.0, 9.0});
    x3bar = x3.serialize();

    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 6.0);
    CHECK(x3bar(1) == 7.0);
    CHECK(x3bar(2) == 8.0);
}

TEST_CASE("RN get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    RN xzero = RN::from_shape(0);
    xzero.unserialize({});
    Eigen::MatrixXd xzerohat = xzero.get_matrix();

    CHECK(xzerohat.rows() == 0);
    CHECK(xzerohat.cols() == 0);

    RN x0 = RN(0);
    x0.unserialize({});
    Eigen::MatrixXd x0hat = x0.get_matrix();

    REQUIRE(x0hat.rows() == 1);
    REQUIRE(x0hat.cols() == 1);
    CHECK(x0hat(0, 0) == 1.0);

    RN x1 = RN(1);
    x1.unserialize({1.0, 2.0});
    Eigen::MatrixXd x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 1) == 1.0);

    CHECK(x1hat(0, 0) == 1.0);
    CHECK(x1hat(1, 0) == 0.0);
    CHECK(x1hat(1, 1) == 1.0);

    x1.unserialize({3.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 1) == 3.0);

    CHECK(x1hat(0, 0) == 1.0);
    CHECK(x1hat(1, 0) == 0.0);
    CHECK(x1hat(1, 1) == 1.0);

    x1.unserialize({});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 1) == 3.0);

    CHECK(x1hat(0, 0) == 1.0);
    CHECK(x1hat(1, 0) == 0.0);
    CHECK(x1hat(1, 1) == 1.0);

    RN x2 = RN(2);
    x2.unserialize({1.0});
    Eigen::MatrixXd x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 2) == 1.0);
    CHECK(x2hat(1, 2) == 0.0);

    CHECK(x2hat(0, 0) == 1.0);
    CHECK(x2hat(0, 1) == 0.0);
    CHECK(x2hat(1, 0) == 0.0);
    CHECK(x2hat(1, 1) == 1.0);
    CHECK(x2hat(2, 0) == 0.0);
    CHECK(x2hat(2, 1) == 0.0);
    CHECK(x2hat(2, 2) == 1.0);

    x2.unserialize({2.0, 3.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 2) == 2.0);
    CHECK(x2hat(1, 2) == 3.0);

    CHECK(x2hat(0, 0) == 1.0);
    CHECK(x2hat(0, 1) == 0.0);
    CHECK(x2hat(1, 0) == 0.0);
    CHECK(x2hat(1, 1) == 1.0);
    CHECK(x2hat(2, 0) == 0.0);
    CHECK(x2hat(2, 1) == 0.0);
    CHECK(x2hat(2, 2) == 1.0);

    x2.unserialize({4.0, 5.0, 6.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 2) == 4.0);
    CHECK(x2hat(1, 2) == 5.0);

    CHECK(x2hat(0, 0) == 1.0);
    CHECK(x2hat(0, 1) == 0.0);
    CHECK(x2hat(1, 0) == 0.0);
    CHECK(x2hat(1, 1) == 1.0);
    CHECK(x2hat(2, 0) == 0.0);
    CHECK(x2hat(2, 1) == 0.0);
    CHECK(x2hat(2, 2) == 1.0);
}

TEST_CASE("RN operator()", "[domain]")
{
    using namespace Lielab::domain;

    RN xzero = RN::from_shape(0);
    xzero.unserialize({});

    // Out of bounds
    CHECK(std::isnan(xzero(-1)));
    CHECK(std::isnan(xzero(0)));
    CHECK(std::isnan(xzero(1)));

    // Out of bounds
    CHECK(std::isnan(xzero(0, -1)));
    CHECK(std::isnan(xzero(-1, 0)));
    CHECK(std::isnan(xzero(-1, -1)));
    CHECK(std::isnan(xzero(0, 0)));
    CHECK(std::isnan(xzero(0, 1)));
    CHECK(std::isnan(xzero(1, 0)));
    CHECK(std::isnan(xzero(1, 1)));

    RN x1 = RN(1);
    x1.unserialize({1.0});

    // In bounds
    CHECK(x1(0) == 1.0);
    CHECK(x1(-1) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x1(-2)));
    CHECK(std::isnan(x1(1)));

    // In bounds
    CHECK(x1(0, 0) == 1.0);
    CHECK(x1(0, 1) == 1.0);
    CHECK(x1(1, 0) == 0.0);
    CHECK(x1(1, 1) == 1.0);
    CHECK(x1(-1, -1) == 1.0);
    CHECK(x1(-1, -2) == 0.0);
    CHECK(x1(-2, -1) == 1.0);
    CHECK(x1(-2, -2) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x1(0, -3)));
    CHECK(std::isnan(x1(-3, 0)));
    CHECK(std::isnan(x1(-3, -3)));
    CHECK(std::isnan(x1(0, 2)));
    CHECK(std::isnan(x1(2, 0)));
    CHECK(std::isnan(x1(2, 2)));

    RN x2 = RN(2);
    x2.unserialize({1.0, 2.0});

    // In bounds
    CHECK(x2(0) == 1.0);
    CHECK(x2(1) == 2.0);
    CHECK(x2(-1) == 2.0);
    CHECK(x2(-2) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x2(-3)));
    CHECK(std::isnan(x2(2)));

    // In bounds
    CHECK(x2(0, 0) == 1.0);
    CHECK(x2(0, 1) == 0.0);
    CHECK(x2(0, 2) == 1.0);
    CHECK(x2(1, 0) == 0.0);
    CHECK(x2(1, 1) == 1.0);
    CHECK(x2(1, 2) == 2.0);
    CHECK(x2(2, 0) == 0.0);
    CHECK(x2(2, 1) == 0.0);
    CHECK(x2(2, 2) == 1.0);
    CHECK(x2(-1, -1) == 1.0);
    CHECK(x2(-1, -2) == 0.0);
    CHECK(x2(-1, -3) == 0.0);
    CHECK(x2(-2, -1) == 2.0);
    CHECK(x2(-2, -2) == 1.0);
    CHECK(x2(-2, -3) == 0.0);
    CHECK(x2(-3, -1) == 1.0);
    CHECK(x2(-3, -2) == 0.0);
    CHECK(x2(-3, -3) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x2(0, -4)));
    CHECK(std::isnan(x2(-4, 0)));
    CHECK(std::isnan(x2(-4, -4)));
    CHECK(std::isnan(x2(0, 3)));
    CHECK(std::isnan(x2(3, 0)));
    CHECK(std::isnan(x2(3, 3)));
}

TEST_CASE("RN math_ops_RN", "[domain]")
{
    using namespace Lielab::domain;

    RN x1(4), x2(4);
    x1.unserialize({1.0, 2.0, 3.0, 4.0});
    x2.unserialize({1.25, 2.5, 3.75, 1.0});

    const RN x1_prod_x2 = x1*x2;
    const Eigen::VectorXd x1_prod_x2bar = x1_prod_x2.serialize();
    REQUIRE(x1_prod_x2bar.size() == 4);
    CHECK(x1_prod_x2bar(0) == 2.25);
    CHECK(x1_prod_x2bar(1) == 4.5);
    CHECK(x1_prod_x2bar(2) == 6.75);
    CHECK(x1_prod_x2bar(3) == 5.0);

    x1 *= x2;
    const Eigen::VectorXd x1bar = x1.serialize();
    REQUIRE(x1bar.size() == 4);
    CHECK(x1bar(0) == 2.25);
    CHECK(x1bar(1) == 4.5);
    CHECK(x1bar(2) == 6.75);
    CHECK(x1bar(3) == 5.0);

    x1.unserialize({1.0, 2.0, 3.0, 4.0});

    const RN x1_inv = x1.inverse();
    const Eigen::VectorXd x1_invbar = x1_inv.serialize();
    REQUIRE(x1_invbar.size() == 4);
    CHECK(x1_invbar(0) == -1.0);
    CHECK(x1_invbar(1) == -2.0);
    CHECK(x1_invbar(2) == -3.0);
    CHECK(x1_invbar(3) == -4.0);
}

TEST_CASE("RN from_vector", "[domain]")
{
    /*!
    * Tests the from_vector operation.
    */

    using namespace Lielab::domain;

    const RN x0 = RN::from_vector({});
    const Eigen::VectorXd x0bar = x0.serialize();

    CHECK(x0.get_shape() == 1);
    CHECK(x0bar.size() == 0);

    const RN x1 = RN::from_vector({1.0});
    const Eigen::VectorXd x1bar = x1.serialize();

    CHECK(x1.get_shape() == 2);
    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == 1.0);

    const RN x2 = RN::from_vector({1.0, 2.0});
    const Eigen::VectorXd x2bar = x2.serialize();

    CHECK(x2.get_shape() == 3);
    REQUIRE(x2bar.size() == 2);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);

    const RN x3 = RN::from_vector({1.0, 2.0, 3.0});
    const Eigen::VectorXd x3bar = x3.serialize();

    CHECK(x3.get_shape() == 4);
    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
}

TEST_CASE("RN project", "[domain]")
{
    using namespace Lielab::domain;

    const Eigen::MatrixXd rand_2_2 = Eigen::MatrixXd::Random(2, 2);
    const Eigen::MatrixXd proj_2_2 = RN::project(rand_2_2);

    REQUIRE(proj_2_2.rows() == 2);
    REQUIRE(proj_2_2.cols() == 2);
    CHECK(proj_2_2(0, 0) == 1.0);
    CHECK(proj_2_2(0, 1) == rand_2_2(0, 1));
    CHECK(proj_2_2(1, 0) == 0.0);
    CHECK(proj_2_2(1, 1) == 1.0);

    const Eigen::MatrixXd rand_3_3 = Eigen::MatrixXd::Random(3, 3);
    const Eigen::MatrixXd proj_3_3 = RN::project(rand_3_3);

    REQUIRE(proj_3_3.rows() == 3);
    REQUIRE(proj_3_3.cols() == 3);
    CHECK(proj_3_3(0, 0) == 1.0);
    CHECK(proj_3_3(0, 1) == 0.0);
    CHECK(proj_3_3(0, 2) == rand_3_3(0, 2));
    CHECK(proj_3_3(1, 0) == 0.0);
    CHECK(proj_3_3(1, 1) == 1.0);
    CHECK(proj_3_3(1, 2) == rand_3_3(1, 2));
    CHECK(proj_3_3(2, 0) == 0.0);
    CHECK(proj_3_3(2, 1) == 0.0);
    CHECK(proj_3_3(2, 2) == 1.0);

    const Eigen::MatrixXd rand_2_3 = Eigen::MatrixXd::Random(2, 3);
    const Eigen::MatrixXd proj_2_3 = RN::project(rand_2_3);

    REQUIRE(proj_2_3.rows() == 2);
    REQUIRE(proj_2_3.cols() == 2);
    CHECK(proj_2_3(0, 0) == 1.0);
    CHECK(proj_2_3(0, 1) == rand_2_3(0, 1));
    CHECK(proj_2_3(1, 0) == 0.0);
    CHECK(proj_2_3(1, 1) == 1.0);

    const Eigen::MatrixXd rand_3_2 = Eigen::MatrixXd::Random(3, 2);
    const Eigen::MatrixXd proj_3_2 = RN::project(rand_3_2);

    REQUIRE(proj_3_2.rows() == 2);
    REQUIRE(proj_3_2.cols() == 2);
    CHECK(proj_3_2(0, 0) == 1.0);
    CHECK(proj_3_2(0, 1) == rand_3_2(0, 1));
    CHECK(proj_3_2(1, 0) == 0.0);
    CHECK(proj_3_2(1, 1) == 1.0);
}
