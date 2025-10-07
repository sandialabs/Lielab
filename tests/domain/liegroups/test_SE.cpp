#include <cmath>
#include <numbers>
#include <Lielab.hpp>
#include <catch2/catch_all.hpp>

TEST_CASE("SE to_string", "[domain]")
{
    using namespace Lielab::domain;

    const SE x0 = SE(0);
    CHECK(x0.to_string() == "SE(0)");
    const SE x1 = SE(1);
    CHECK(x1.to_string() == "SE(1)");
    const SE x10 = SE(10);
    CHECK(x10.to_string() == "SE(10)");
}

TEST_CASE("SE main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const SE xblank = SE();
    CHECK(xblank.get_dimension() == 0);

    const SE x0 = SE(0);
    CHECK(x0.get_dimension() == 0);
    const SE x1 = SE(1);
    CHECK(x1.get_dimension() == 1);
    const SE x10 = SE(10);
    CHECK(x10.get_dimension() == 55);
}

TEST_CASE("SE matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const SE x0 = SE(Eigen::MatrixXd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    const SE x1 = SE(Eigen::MatrixXd::Random(1, 1));
    CHECK(x1.get_shape() == 1);

    const SE x2 = SE(Eigen::MatrixXd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(SE(Eigen::MatrixXd::Random(2, 3)));
    CHECK_THROWS(SE(Eigen::MatrixXd::Random(3, 2)));
}

TEST_CASE("SE from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const SE x0 = SE::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const SE x1 = SE::from_shape(1);
    CHECK(x1.get_dimension() == 0);
    const Eigen::MatrixXd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const SE x2 = SE::from_shape(2);
    CHECK(x2.get_dimension() == 1);
    const Eigen::MatrixXd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("SE get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    SE veryzero = SE::from_shape(0);
    SE zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(veryzero.get_dimension() == 0);
    CHECK(zero.get_dimension() == 0);
    CHECK(one.get_dimension() == 1);
    CHECK(two.get_dimension() == 3);
    CHECK(three.get_dimension() == 6);
    CHECK(four.get_dimension() == 10);
    CHECK(five.get_dimension() == 15);
    CHECK(six.get_dimension() == 21);
    CHECK(seven.get_dimension() == 28);
    CHECK(eight.get_dimension() == 36);
}

TEST_CASE("SE serialize/unserialize", "[domain]")
{
    /*!
    * Tests the serialize/unserialize operation.
    */

    using namespace Lielab::domain;

    SE xzero = SE::from_shape(0);
    xzero.unserialize({});
    Eigen::VectorXd xzerobar = xzero.serialize();

    CHECK(xzerobar.size() == 0);

    SE x0 = SE(0);
    x0.unserialize({});
    Eigen::VectorXd x0bar = x0.serialize();

    CHECK(x0bar.size() == 1);
    CHECK(x0bar(0) == 1.0);

    SE x2 = SE(2);
    x2.unserialize({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});
    Eigen::VectorXd x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 9);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 3.0);
    CHECK(x2bar(3) == 4.0);
    CHECK(x2bar(4) == 5.0);
    CHECK(x2bar(5) == 6.0);
    CHECK(x2bar(6) == 7.0);
    CHECK(x2bar(7) == 8.0);
    CHECK(x2bar(8) == 9.0);

    x2.unserialize({10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0});
    x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 9);
    CHECK(x2bar(0) == 10.0);
    CHECK(x2bar(1) == 11.0);
    CHECK(x2bar(2) == 12.0);
    CHECK(x2bar(3) == 13.0);
    CHECK(x2bar(4) == 14.0);
    CHECK(x2bar(5) == 15.0);
    CHECK(x2bar(6) == 16.0);
    CHECK(x2bar(7) == 17.0);
    CHECK(x2bar(8) == 18.0);

    x2.unserialize({19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0});
    x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 9);
    CHECK(x2bar(0) == 19.0);
    CHECK(x2bar(1) == 20.0);
    CHECK(x2bar(2) == 21.0);
    CHECK(x2bar(3) == 22.0);
    CHECK(x2bar(4) == 23.0);
    CHECK(x2bar(5) == 24.0);
    CHECK(x2bar(6) == 25.0);
    CHECK(x2bar(7) == 26.0);
    CHECK(x2bar(8) == 18.0);

    SE x3 = SE(3);
    x3.unserialize({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0});
    Eigen::VectorXd x3bar = x3.serialize();

    REQUIRE(x3bar.size() == 16);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
    CHECK(x3bar(3) == 4.0);
    CHECK(x3bar(4) == 5.0);
    CHECK(x3bar(5) == 6.0);
    CHECK(x3bar(6) == 7.0);
    CHECK(x3bar(7) == 8.0);
    CHECK(x3bar(8) == 9.0);
    CHECK(x3bar(9) == 10.0);
    CHECK(x3bar(10) == 11.0);
    CHECK(x3bar(11) == 12.0);
    CHECK(x3bar(12) == 13.0);
    CHECK(x3bar(13) == 14.0);
    CHECK(x3bar(14) == 15.0);
    CHECK(x3bar(15) == 1.0);

    x3.unserialize({16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0});
    x3bar = x3.serialize();

    REQUIRE(x3bar.size() == 16);
    CHECK(x3bar(0) == 16.0);
    CHECK(x3bar(1) == 17.0);
    CHECK(x3bar(2) == 18.0);
    CHECK(x3bar(3) == 19.0);
    CHECK(x3bar(4) == 20.0);
    CHECK(x3bar(5) == 21.0);
    CHECK(x3bar(6) == 22.0);
    CHECK(x3bar(7) == 23.0);
    CHECK(x3bar(8) == 24.0);
    CHECK(x3bar(9) == 25.0);
    CHECK(x3bar(10) == 26.0);
    CHECK(x3bar(11) == 27.0);
    CHECK(x3bar(12) == 28.0);
    CHECK(x3bar(13) == 29.0);
    CHECK(x3bar(14) == 30.0);
    CHECK(x3bar(15) == 31.0);

    x3.unserialize({32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0});
    x3bar = x3.serialize();

    REQUIRE(x3bar.size() == 16);
    CHECK(x3bar(0) == 32.0);
    CHECK(x3bar(1) == 33.0);
    CHECK(x3bar(2) == 34.0);
    CHECK(x3bar(3) == 35.0);
    CHECK(x3bar(4) == 36.0);
    CHECK(x3bar(5) == 37.0);
    CHECK(x3bar(6) == 38.0);
    CHECK(x3bar(7) == 39.0);
    CHECK(x3bar(8) == 40.0);
    CHECK(x3bar(9) == 41.0);
    CHECK(x3bar(10) == 42.0);
    CHECK(x3bar(11) == 43.0);
    CHECK(x3bar(12) == 44.0);
    CHECK(x3bar(13) == 45.0);
    CHECK(x3bar(14) == 46.0);
    CHECK(x3bar(15) == 47.0);
}

TEST_CASE("SE get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    SE xzero = SE::from_shape(0);
    xzero.unserialize({});
    Eigen::MatrixXd xzerohat = xzero.get_matrix();

    CHECK(xzerohat.rows() == 0);
    CHECK(xzerohat.cols() == 0);

    SE x0 = SE(0);
    x0.unserialize({});
    Eigen::MatrixXd x0hat = x0.get_matrix();

    REQUIRE(x0hat.rows() == 1);
    REQUIRE(x0hat.cols() == 1);
    CHECK(x0hat(0, 0) == 1.0);

    SE x1 = SE(1);
    x1.unserialize({1.0, 2.0});
    Eigen::MatrixXd x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 0) == 1.0);
    CHECK(x1hat(0, 1) == 2.0);
    CHECK(x1hat(1, 0) == 0.0);
    CHECK(x1hat(1, 1) == 1.0);

    x1.unserialize({3.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 0) == 3.0);
    CHECK(x1hat(0, 1) == 2.0);
    CHECK(x1hat(1, 0) == 0.0);
    CHECK(x1hat(1, 1) == 1.0);

    x1.unserialize({});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 0) == 3.0);
    CHECK(x1hat(0, 1) == 2.0);
    CHECK(x1hat(1, 0) == 0.0);
    CHECK(x1hat(1, 1) == 1.0);

    SE x2 = SE(2);
    x2.unserialize({1.0});
    Eigen::MatrixXd x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 0) == 1.0);
    CHECK(x2hat(0, 1) == 0.0);
    CHECK(x2hat(0, 2) == 0.0);
    CHECK(x2hat(1, 0) == 0.0);
    CHECK(x2hat(1, 1) == 1.0);
    CHECK(x2hat(1, 2) == 0.0);
    CHECK(x2hat(2, 0) == 0.0);
    CHECK(x2hat(2, 1) == 0.0);
    CHECK(x2hat(2, 2) == 1.0);

    x2.unserialize({2.0, 3.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 0) == 2.0);
    CHECK(x2hat(0, 1) == 3.0);
    CHECK(x2hat(0, 2) == 0.0);
    CHECK(x2hat(1, 0) == 0.0);
    CHECK(x2hat(1, 1) == 1.0);
    CHECK(x2hat(1, 2) == 0.0);
    CHECK(x2hat(2, 0) == 0.0);
    CHECK(x2hat(2, 1) == 0.0);
    CHECK(x2hat(2, 2) == 1.0);

    x2.unserialize({4.0, 5.0, 6.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 0) == 4.0);
    CHECK(x2hat(0, 1) == 5.0);
    CHECK(x2hat(0, 2) == 6.0);
    CHECK(x2hat(1, 0) == 0.0);
    CHECK(x2hat(1, 1) == 1.0);
    CHECK(x2hat(1, 2) == 0.0);
    CHECK(x2hat(2, 0) == 0.0);
    CHECK(x2hat(2, 1) == 0.0);
    CHECK(x2hat(2, 2) == 1.0);
}

TEST_CASE("SE operator()", "[domain]")
{
    using namespace Lielab::domain;

    SE xzero = SE::from_shape(0);
    xzero.unserialize({});

    // Out of bounds
    // CHECK(std::isnan(xzero(-1)));
    // CHECK(std::isnan(xzero(0)));
    // CHECK(std::isnan(xzero(1)));

    // Out of bounds
    CHECK(std::isnan(xzero(0, -1)));
    CHECK(std::isnan(xzero(-1, 0)));
    CHECK(std::isnan(xzero(-1, -1)));
    CHECK(std::isnan(xzero(0, 0)));
    CHECK(std::isnan(xzero(0, 1)));
    CHECK(std::isnan(xzero(1, 0)));
    CHECK(std::isnan(xzero(1, 1)));

    SE x1 = SE(1);
    x1.unserialize({1.0});

    // In bounds
    // CHECK(x1(0) == 1.0);

    // Out of bounds
    // CHECK(std::isnan(x1(-1)));
    // CHECK(std::isnan(x1(1)));

    // In bounds
    CHECK(x1(0, 0) == 1.0);
    CHECK(x1(0, 1) == 0.0);
    CHECK(x1(1, 0) == 0.0);
    CHECK(x1(1, 1) == 1.0);
    CHECK(x1(-1, -1) == 1.0);
    CHECK(x1(-1, -2) == 0.0);
    CHECK(x1(-2, -1) == 0.0);
    CHECK(x1(-2, -2) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x1(0, -3)));
    CHECK(std::isnan(x1(-3, 0)));
    CHECK(std::isnan(x1(-3, -3)));
    CHECK(std::isnan(x1(0, 2)));
    CHECK(std::isnan(x1(2, 0)));
    CHECK(std::isnan(x1(2, 2)));

    SE x2 = SE(2);
    x2.unserialize({1.0, 2.0, 3.0});

    // In bounds
    // CHECK(x2(0) == 1.0);
    // CHECK(x2(1) == 2.0);
    // CHECK(x2(2) == 3.0);

    // Out of bounds
    // CHECK(std::isnan(x2(-1)));
    // CHECK(std::isnan(x2(3)));

    // In bounds
    CHECK(x2(0, 0) == 1.0);
    CHECK(x2(0, 1) == 2.0);
    CHECK(x2(0, 2) == 3.0);
    CHECK(x2(1, 0) == 0.0);
    CHECK(x2(1, 1) == 1.0);
    CHECK(x2(1, 2) == 0.0);
    CHECK(x2(2, 0) == 0.0);
    CHECK(x2(2, 1) == 0.0);
    CHECK(x2(2, 2) == 1.0);
    CHECK(x2(-1, -1) == 1.0);
    CHECK(x2(-1, -2) == 0.0);
    CHECK(x2(-1, -3) == 0.0);
    CHECK(x2(-2, -1) == 0.0);
    CHECK(x2(-2, -2) == 1.0);
    CHECK(x2(-2, -3) == 0.0);
    CHECK(x2(-3, -1) == 3.0);
    CHECK(x2(-3, -2) == 2.0);
    CHECK(x2(-3, -3) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x2(0, -4)));
    CHECK(std::isnan(x2(-4, 0)));
    CHECK(std::isnan(x2(-4, -4)));
    CHECK(std::isnan(x2(0, 3)));
    CHECK(std::isnan(x2(3, 0)));
    CHECK(std::isnan(x2(3, 3)));
}

TEST_CASE("SE math_ops_SE", "[domain]")
{
    using namespace Lielab::domain;

    SE x1(2), x2(2);
    x1.unserialize({1.0, 2.0, 3.0});
    x2.unserialize({1.25, 2.5, 3.75});

    const SE x1_prod_x2 = x1*x2;
    const Eigen::VectorXd x1_prod_x2bar = x1_prod_x2.serialize();
    REQUIRE(x1_prod_x2bar.size() == 9);
    CHECK(x1_prod_x2bar(0) == 1.25);
    CHECK(x1_prod_x2bar(1) == 4.5);
    CHECK(x1_prod_x2bar(2) == 6.75);
    CHECK(x1_prod_x2bar(3) == 0.0);
    CHECK(x1_prod_x2bar(4) == 1.0);
    CHECK(x1_prod_x2bar(5) == 0.0);
    CHECK(x1_prod_x2bar(6) == 0.0);
    CHECK(x1_prod_x2bar(7) == 0.0);
    CHECK(x1_prod_x2bar(8) == 1.0);

    x1 *= x2;
    const Eigen::VectorXd x1bar = x1.serialize();
    REQUIRE(x1bar.size() == 9);
    CHECK(x1bar(0) == 1.25);
    CHECK(x1bar(1) == 4.5);
    CHECK(x1bar(2) == 6.75);
    CHECK(x1bar(3) == 0.0);
    CHECK(x1bar(4) == 1.0);
    CHECK(x1bar(5) == 0.0);
    CHECK(x1bar(6) == 0.0);
    CHECK(x1bar(7) == 0.0);
    CHECK(x1bar(8) == 1.0);

    x1.unserialize({1.0, 2.0, 3.0});

    const SE x1_inv = x1.inverse();
    const Eigen::VectorXd x1_invbar = x1_inv.serialize();
    REQUIRE(x1_invbar.size() == 9);
    CHECK_THAT(x1_invbar(0), Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(x1_invbar(1), Catch::Matchers::WithinULP(-2.0, 0));
    CHECK_THAT(x1_invbar(2), Catch::Matchers::WithinULP(-3.0, 0));
    CHECK_THAT(x1_invbar(3), Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(x1_invbar(4), Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(x1_invbar(5), Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(x1_invbar(6), Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(x1_invbar(7), Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(x1_invbar(8), Catch::Matchers::WithinULP(1.0, 0));
}

// TEST_CASE("SE project", "[domain]")
// {
//     using namespace Lielab::domain;

//     const Eigen::MatrixXd rand_2_2 = Eigen::MatrixXd::Random(2, 2);
//     const Eigen::MatrixXd proj_2_2 = SO::project(rand_2_2);

//     REQUIRE(proj_2_2.rows() == 2);
//     REQUIRE(proj_2_2.cols() == 2);
//     CHECK_THAT(std::abs(proj_2_2.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));

//     const Eigen::MatrixXd rand_3_3 = Eigen::MatrixXd::Random(3, 3);
//     const Eigen::MatrixXd proj_3_3 = SO::project(rand_3_3);

//     REQUIRE(proj_3_3.rows() == 3);
//     REQUIRE(proj_3_3.cols() == 3);
//     CHECK_THAT(std::abs(proj_3_3.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));

//     const Eigen::MatrixXd rand_2_3 = Eigen::MatrixXd::Random(2, 3);
//     const Eigen::MatrixXd proj_2_3 = SO::project(rand_2_3);

//     REQUIRE(proj_2_3.rows() == 2);
//     REQUIRE(proj_2_3.cols() == 2);
//     CHECK_THAT(std::abs(proj_2_3.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));

//     const Eigen::MatrixXd rand_3_2 = Eigen::MatrixXd::Random(3, 2);
//     const Eigen::MatrixXd proj_3_2 = SO::project(rand_3_2);

//     REQUIRE(proj_3_2.rows() == 2);
//     REQUIRE(proj_3_2.cols() == 2);
//     CHECK_THAT(std::abs(proj_3_2.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));
// }

