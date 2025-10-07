#include <cmath>
#include <numbers>
#include <Lielab.hpp>
#include <catch2/catch_all.hpp>

TEST_CASE("SO to_string", "[domain]")
{
    using namespace Lielab::domain;

    const SO x0 = SO(0);
    CHECK(x0.to_string() == "SO(0)");
    const SO x1 = SO(1);
    CHECK(x1.to_string() == "SO(1)");
    const SO x10 = SO(10);
    CHECK(x10.to_string() == "SO(10)");
}

TEST_CASE("SO main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const SO xblank = SO();
    CHECK(xblank.get_dimension() == 0);

    const SO x0 = SO(0);
    CHECK(x0.get_dimension() == 0);
    const SO x1 = SO(1);
    CHECK(x1.get_dimension() == 0);
    const SO x10 = SO(10);
    CHECK(x10.get_dimension() == 45);
}

TEST_CASE("SO matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const SO x0 = SO(Eigen::MatrixXd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    const SO x1 = SO(Eigen::MatrixXd::Random(1, 1));
    CHECK(x1.get_shape() == 1);

    const SO x2 = SO(Eigen::MatrixXd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(SO(Eigen::MatrixXd::Random(2, 3)));
    CHECK_THROWS(SO(Eigen::MatrixXd::Random(3, 2)));
}

TEST_CASE("SO from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const SO x0 = SO::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const SO x1 = SO::from_shape(1);
    CHECK(x1.get_dimension() == 0);
    const Eigen::MatrixXd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const SO x2 = SO::from_shape(2);
    CHECK(x2.get_dimension() == 1);
    const Eigen::MatrixXd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("SO get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    SO zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

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

TEST_CASE("SO serialize/unserialize", "[domain]")
{
    /*!
    * Tests the serialize/unserialize operation.
    */

    using namespace Lielab::domain;

    SO x0 = SO(0);
    x0.unserialize({});
    Eigen::VectorXd x0bar = x0.serialize();

    CHECK(x0bar.size() == 0);

    SO x1 = SO(1);
    x1.unserialize({1.0, 2.0});
    Eigen::VectorXd x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == 1.0);

    x1.unserialize({3.0});
    x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == 3.0);

    x1.unserialize({});
    x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == 3.0);

    SO x2 = SO(2);
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

TEST_CASE("SO get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    SO x0 = SO(0);
    x0.unserialize({});
    Eigen::MatrixXd x0hat = x0.get_matrix();

    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    SO x1 = SO(1);
    x1.unserialize({});
    Eigen::MatrixXd x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == 1.0);

    x1.unserialize({1.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == 1.0);

    x1.unserialize({2.0, 3.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == 2.0);

    SO x2 = SO(2);
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

TEST_CASE("SO operator()", "[domain]")
{
    using namespace Lielab::domain;

    SO x0 = SO::from_shape(0);
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

    SO x1 = SO(1);
    x1.unserialize({1.0});

    // In bounds
    // CHECK(x1(0) == 1.0);

    // Out of bounds
    // CHECK(std::isnan(x1(-1)));
    // CHECK(std::isnan(x1(1)));

    // In bounds
    CHECK(x1(0, 0) == 1.0);
    CHECK(x1(-1, -1) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x1(0, -2)));
    CHECK(std::isnan(x1(-2, 0)));
    CHECK(std::isnan(x1(-2, -2)));
    CHECK(std::isnan(x1(0, 1)));
    CHECK(std::isnan(x1(1, 0)));
    CHECK(std::isnan(x1(1, 1)));

    SO x2 = SO(2);
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

TEST_CASE("SO math_ops_SO", "[domain]")
{
    using namespace Lielab::domain;

    SO x1(2), x2(2);
    x1.unserialize({1.0, 2.0, 3.0, 4.0});
    x2.unserialize({1.25, 2.5, 3.75, 1.0});

    const SO x1_prod_x2 = x1*x2;
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

    const SO x1_inv = x1.inverse();
    const Eigen::VectorXd x1_invbar = x1_inv.serialize();
    REQUIRE(x1_invbar.size() == 4);
    CHECK_THAT(x1_invbar(0), Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(x1_invbar(1), Catch::Matchers::WithinULP(3.0, 0));
    CHECK_THAT(x1_invbar(2), Catch::Matchers::WithinULP(2.0, 0));
    CHECK_THAT(x1_invbar(3), Catch::Matchers::WithinULP(4.0, 0));
}

TEST_CASE("SO project", "[domain]")
{
    using namespace Lielab::domain;

    const Eigen::MatrixXd rand_2_2 = Eigen::MatrixXd::Random(2, 2);
    const Eigen::MatrixXd proj_2_2 = SO::project(rand_2_2);

    REQUIRE(proj_2_2.rows() == 2);
    REQUIRE(proj_2_2.cols() == 2);
    CHECK_THAT(std::abs(proj_2_2.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));

    const Eigen::MatrixXd rand_3_3 = Eigen::MatrixXd::Random(3, 3);
    const Eigen::MatrixXd proj_3_3 = SO::project(rand_3_3);

    REQUIRE(proj_3_3.rows() == 3);
    REQUIRE(proj_3_3.cols() == 3);
    CHECK_THAT(std::abs(proj_3_3.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));

    const Eigen::MatrixXd rand_2_3 = Eigen::MatrixXd::Random(2, 3);
    const Eigen::MatrixXd proj_2_3 = SO::project(rand_2_3);

    REQUIRE(proj_2_3.rows() == 2);
    REQUIRE(proj_2_3.cols() == 2);
    CHECK_THAT(std::abs(proj_2_3.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));

    const Eigen::MatrixXd rand_3_2 = Eigen::MatrixXd::Random(3, 2);
    const Eigen::MatrixXd proj_3_2 = SO::project(rand_3_2);

    REQUIRE(proj_3_2.rows() == 2);
    REQUIRE(proj_3_2.cols() == 2);
    CHECK_THAT(std::abs(proj_3_2.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14));
}

