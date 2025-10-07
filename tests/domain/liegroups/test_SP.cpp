#include <Lielab.hpp>
#include <iostream>

#include <catch2/catch_all.hpp>

TEST_CASE("SP to_string", "[domain]")
{
    using namespace Lielab::domain;

    const SP x0 = SP(0);
    CHECK(x0.to_string() == "SP(0, R)");
    CHECK_THROWS(SP(1));
    const SP x10 = SP(10);
    CHECK(x10.to_string() == "SP(10, R)");
}

TEST_CASE("SP main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const SP xblank = SP();
    CHECK(xblank.get_dimension() == 0);

    const SP x0 = SP(0);
    CHECK(x0.get_dimension() == 0);
    CHECK_THROWS(SP(1));
    const SP x10 = SP(10);
    CHECK(x10.get_dimension() == 55);
}

TEST_CASE("SP matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const SP x0 = SP(Eigen::MatrixXd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    CHECK_THROWS(SP(Eigen::MatrixXd::Random(1, 1)));

    const SP x2 = SP(Eigen::MatrixXd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(SP(Eigen::MatrixXd::Random(2, 3)));
    CHECK_THROWS(SP(Eigen::MatrixXd::Random(3, 2)));
}

TEST_CASE("SP from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const SP x0 = SP::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    CHECK_THROWS(SP::from_shape(1));

    const SP x2 = SP::from_shape(2);
    CHECK(x2.get_dimension() == 3);
    const Eigen::MatrixXd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("SP get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    SP zero(0), two(2), four(4), six(6), eight(8);

    CHECK(zero.get_dimension() == 0);
    CHECK(two.get_dimension() == 3);
    CHECK(four.get_dimension() == 10);
    CHECK(six.get_dimension() == 21);
    CHECK(eight.get_dimension() == 36);
}

TEST_CASE("SP serialize/unserialize", "[domain]")
{
    /*!
    * Tests the serialize/unserialize operation.
    */

    using namespace Lielab::domain;

    SP x0 = SP(0);
    x0.unserialize({});
    Eigen::VectorXd x0bar = x0.serialize();

    CHECK(x0bar.size() == 0);

    SP x1 = SP(2);
    x1.unserialize({1.0, 2.0, 3.0, 4.0, 5.0});
    Eigen::VectorXd x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 4);
    CHECK(x1bar(0) == 1.0);
    CHECK(x1bar(1) == 2.0);
    CHECK(x1bar(2) == 3.0);
    CHECK(x1bar(3) == 4.0);

    x1.unserialize({6.0, 7.0, 8.0, 9.0});
    x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 4);
    CHECK(x1bar(0) == 6.0);
    CHECK(x1bar(1) == 7.0);
    CHECK(x1bar(2) == 8.0);
    CHECK(x1bar(3) == 9.0);

    x1.unserialize({10.0, 11.0, 12.0});
    x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 4);
    CHECK(x1bar(0) == 10.0);
    CHECK(x1bar(1) == 11.0);
    CHECK(x1bar(2) == 12.0);
    CHECK(x1bar(3) == 9.0);

    SP x2 = SP(2);
    x2.unserialize({1.0, 2.0, 3.0});
    Eigen::VectorXd x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 3.0);
    CHECK(x2bar(3) == 1.0);

    x2.unserialize({4.0, 5.0, 6.0, 7.0});
    x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 4.0);
    CHECK(x2bar(1) == 5.0);
    CHECK(x2bar(2) == 6.0);
    CHECK(x2bar(3) == 7.0);

    x2.unserialize({8.0, 9.0, 10.0, 11.0, 12.0});
    x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 8.0);
    CHECK(x2bar(1) == 9.0);
    CHECK(x2bar(2) == 10.0);
    CHECK(x2bar(3) == 11.0);
}

TEST_CASE("SP get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    SP x0 = SP(0);
    x0.unserialize({});
    Eigen::MatrixXd x0hat = x0.get_matrix();

    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    SP x2 = SP(2);
    x2.unserialize({1.0, 2.0, 3.0});
    Eigen::MatrixXd x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 1.0);
    CHECK(x2hat(0, 1) == 2.0);
    CHECK(x2hat(1, 0) == 3.0);
    CHECK(x2hat(1, 1) == 1.0);

    x2.unserialize({4.0, 5.0, 6.0, 7.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 4.0);
    CHECK(x2hat(0, 1) == 5.0);
    CHECK(x2hat(1, 0) == 6.0);
    CHECK(x2hat(1, 1) == 7.0);

    x2.unserialize({8.0, 9.0, 10.0, 11.0, 12.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == 8.0);
    CHECK(x2hat(0, 1) == 9.0);
    CHECK(x2hat(1, 0) == 10.0);
    CHECK(x2hat(1, 1) == 11.0);
}

TEST_CASE("SP operator()", "[domain]")
{
    using namespace Lielab::domain;

    SP x0 = SP::from_shape(0);
    x0.unserialize({});

    // Out of bounds
    // CHECK(std::isnan(x0(-1)));
    // CHECK(std::isnan(x0(0)));
    // CHECK(std::isnan(x0(1)));

    // Out of bounds
    CHECK(std::isnan(x0(0, -1)));
    CHECK(std::isnan(x0(-1, 0)));
    CHECK(std::isnan(x0(-1, -1)));
    CHECK(std::isnan(x0(0, 0)));
    CHECK(std::isnan(x0(0, 1)));
    CHECK(std::isnan(x0(1, 0)));
    CHECK(std::isnan(x0(1, 1)));

    SP x2 = SP(2);
    x2.unserialize({1.0, 2.0, 3.0, 4.0});

    // In bounds
    // CHECK(x2(0) == 1.0);
    // CHECK(x2(1) == 2.0);
    // CHECK(x2(2) == 3.0);
    // CHECK(x2(3) == 4.0);

    // Out of bounds
    // CHECK(std::isnan(x2(-1)));
    // CHECK(std::isnan(x2(4)));

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

TEST_CASE("SP math_ops_SP", "[domain]")
{
    using namespace Lielab::domain;

    SP x1(2), x2(2);
    x1.unserialize({1.0, 2.0, 3.0, 4.0});
    x2.unserialize({1.25, 2.5, 3.75, 1.0});

    const SP x1_prod_x2 = x1*x2;
    const Eigen::VectorXd x1_prod_x2bar = x1_prod_x2.serialize();
    REQUIRE(x1_prod_x2bar.size() == 4);
    CHECK(x1_prod_x2bar(0) == 8.75);
    CHECK(x1_prod_x2bar(1) == 4.5);
    CHECK(x1_prod_x2bar(2) == 18.75);
    CHECK(x1_prod_x2bar(3) == 11.5);

    x1 *= x2;
    const Eigen::VectorXd x1bar = x1.serialize();
    REQUIRE(x1bar.size() == 4);
    CHECK(x1bar(0) == 8.75);
    CHECK(x1bar(1) == 4.5);
    CHECK(x1bar(2) == 18.75);
    CHECK(x1bar(3) == 11.5);

    x1.unserialize({1.0, 2.0, 3.0, 4.0});

    const SP x1_inv = x1.inverse();
    const Eigen::VectorXd x1_invbar = x1_inv.serialize();
    REQUIRE(x1_invbar.size() == 4);
    CHECK_THAT(x1_invbar(0), Catch::Matchers::WithinULP(-2.0, 2));
    CHECK_THAT(x1_invbar(1), Catch::Matchers::WithinULP(1.0, 2));
    CHECK_THAT(x1_invbar(2), Catch::Matchers::WithinULP(1.5, 2));
    CHECK_THAT(x1_invbar(3), Catch::Matchers::WithinULP(-0.5, 2));
}

// TEST_CASE("SP project", "[domain]")
// {
//     using namespace Lielab::domain;

//     const Eigen::MatrixXd rand_2_2 = Eigen::MatrixXd::Random(2, 2);
//     const Eigen::MatrixXd proj_2_2 = SP::project(rand_2_2);

//     REQUIRE(proj_2_2.rows() == 2);
//     REQUIRE(proj_2_2.cols() == 2);
//     CHECK_THAT(std::abs(proj_2_2.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));

//     const Eigen::MatrixXd rand_3_3 = Eigen::MatrixXd::Random(3, 3);
//     const Eigen::MatrixXd proj_3_3 = SP::project(rand_3_3);

//     REQUIRE(proj_3_3.rows() == 3);
//     REQUIRE(proj_3_3.cols() == 3);
//     CHECK_THAT(std::abs(proj_3_3.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));

//     const Eigen::MatrixXd rand_2_3 = Eigen::MatrixXd::Random(2, 3);
//     const Eigen::MatrixXd proj_2_3 = SP::project(rand_2_3);

//     REQUIRE(proj_2_3.rows() == 2);
//     REQUIRE(proj_2_3.cols() == 2);
//     CHECK_THAT(std::abs(proj_2_3.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));

//     const Eigen::MatrixXd rand_3_2 = Eigen::MatrixXd::Random(3, 2);
//     const Eigen::MatrixXd proj_3_2 = SP::project(rand_3_2);

//     REQUIRE(proj_3_2.rows() == 2);
//     REQUIRE(proj_3_2.cols() == 2);
//     CHECK_THAT(std::abs(proj_3_2.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));
// }

