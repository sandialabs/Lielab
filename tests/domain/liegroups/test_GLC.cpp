#include <Lielab.hpp>
#include <iostream>

#include <catch2/catch_all.hpp>

TEST_CASE("GLC to_string", "[domain]")
{
    using namespace Lielab::domain;

    const GLC x0 = GLC(0);
    CHECK(x0.to_string() == "GL(0, C)");
    const GLC x1 = GLC(1);
    CHECK(x1.to_string() == "GL(1, C)");
    const GLC x10 = GLC(10);
    CHECK(x10.to_string() == "GL(10, C)");
}

TEST_CASE("GLC main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const GLC xblank = GLC();
    CHECK(xblank.get_dimension() == 0);

    const GLC x0 = GLC(0);
    CHECK(x0.get_dimension() == 0);
    const GLC x1 = GLC(1);
    CHECK(x1.get_dimension() == 2);
    const GLC x10 = GLC(10);
    CHECK(x10.get_dimension() == 200);
}

TEST_CASE("GLC matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const GLC x0 = GLC(Eigen::MatrixXcd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    const GLC x1 = GLC(Eigen::MatrixXcd::Random(1, 1));
    CHECK(x1.get_shape() == 1);

    const GLC x2 = GLC(Eigen::MatrixXcd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(GLC(Eigen::MatrixXcd::Random(2, 3)));
    CHECK_THROWS(GLC(Eigen::MatrixXcd::Random(3, 2)));
}

TEST_CASE("GLC from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const GLC x0 = GLC::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXcd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const GLC x1 = GLC::from_shape(1);
    CHECK(x1.get_dimension() == 2);
    const Eigen::MatrixXcd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const GLC x2 = GLC::from_shape(2);
    CHECK(x2.get_dimension() == 8);
    const Eigen::MatrixXcd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("GLC get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    GLC zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    CHECK(zero.get_dimension() == 0);
    CHECK(one.get_dimension() == 2);
    CHECK(two.get_dimension() == 8);
    CHECK(three.get_dimension() == 18);
    CHECK(four.get_dimension() == 32);
    CHECK(five.get_dimension() == 50);
    CHECK(six.get_dimension() == 72);
    CHECK(seven.get_dimension() == 98);
    CHECK(eight.get_dimension() == 128);
}

TEST_CASE("GLC serialize/unserialize", "[domain]")
{
    /*!
    * Tests the serialize/unserialize operation.
    */

    using namespace Lielab::domain;

    GLC x0 = GLC(0);
    x0.unserialize({});
    Eigen::VectorXd x0bar = x0.serialize();

    CHECK(x0bar.size() == 0);

    GLC x1 = GLC(1);
    x1.unserialize({1.0, 2.0, 3.0});
    Eigen::VectorXd x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 1.0);
    CHECK(x1bar(1) == 2.0);

    x1.unserialize({4.0, 5.0});
    x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 4.0);
    CHECK(x1bar(1) == 5.0);

    x1.unserialize({6.0});
    x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 6.0);
    CHECK(x1bar(1) == 5.0);

    GLC x2 = GLC(2);
    x2.unserialize({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});
    Eigen::VectorXd x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 8);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 3.0);
    CHECK(x2bar(3) == 4.0);
    CHECK(x2bar(4) == 5.0);
    CHECK(x2bar(5) == 6.0);
    CHECK(x2bar(6) == 7.0);
    CHECK(x2bar(7) == 0.0);

    x2.unserialize({8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0});
    x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 8);
    CHECK(x2bar(0) == 8.0);
    CHECK(x2bar(1) == 9.0);
    CHECK(x2bar(2) == 10.0);
    CHECK(x2bar(3) == 11.0);
    CHECK(x2bar(4) == 12.0);
    CHECK(x2bar(5) == 13.0);
    CHECK(x2bar(6) == 14.0);
    CHECK(x2bar(7) == 15.0);

    x2.unserialize({16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0});
    x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 8);
    CHECK(x2bar(0) == 16.0);
    CHECK(x2bar(1) == 17.0);
    CHECK(x2bar(2) == 18.0);
    CHECK(x2bar(3) == 19.0);
    CHECK(x2bar(4) == 20.0);
    CHECK(x2bar(5) == 21.0);
    CHECK(x2bar(6) == 22.0);
    CHECK(x2bar(7) == 23.0);
}

TEST_CASE("GLC get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    GLC x0 = GLC(0);
    x0.unserialize({});
    Eigen::MatrixXcd x0hat = x0.get_matrix();

    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    GLC x1 = GLC(1);
    x1.unserialize({1.0, 2.0, 3.0});
    Eigen::MatrixXcd x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == std::complex<double>(1.0, 2.0));

    x1.unserialize({4.0, 5.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == std::complex<double>(4.0, 5.0));

    x1.unserialize({6.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == std::complex<double>(6.0, 5.0));

    GLC x2 = GLC(2);
    x2.unserialize({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});
    Eigen::MatrixXcd x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == std::complex<double>(1.0, 2.0));
    CHECK(x2hat(0, 1) == std::complex<double>(3.0, 4.0));
    CHECK(x2hat(1, 0) == std::complex<double>(5.0, 6.0));
    CHECK(x2hat(1, 1) == std::complex<double>(7.0, 0.0));

    x2.unserialize({8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == std::complex<double>(8.0, 9.0));
    CHECK(x2hat(0, 1) == std::complex<double>(10.0, 11.0));
    CHECK(x2hat(1, 0) == std::complex<double>(12.0, 13.0));
    CHECK(x2hat(1, 1) == std::complex<double>(14.0, 15.0));

    x2.unserialize({16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == std::complex<double>(16.0, 17.0));
    CHECK(x2hat(0, 1) == std::complex<double>(18.0, 19.0));
    CHECK(x2hat(1, 0) == std::complex<double>(20.0, 21.0));
    CHECK(x2hat(1, 1) == std::complex<double>(22.0, 23.0));
}

TEST_CASE("GLC operator()", "[domain]")
{
    using namespace Lielab::domain;

    GLC x0 = GLC::from_shape(0);
    x0.unserialize({});

    // Out of bounds
    // CHECK(std::isnan(x0(-1)));
    // CHECK(std::isnan(x0(0)));
    // CHECK(std::isnan(x0(1)));

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

    GLC x1 = GLC(1);
    x1.unserialize({1.0, 2.0});

    // In bounds
    // CHECK(x1(0) == 1.0);
    // CHECK(x1(1) == 2.0);

    // Out of bounds
    // CHECK(std::isnan(x1(-1)));
    // CHECK(std::isnan(x1(2)));

    // In bounds
    CHECK(x1(0, 0) == std::complex<double>(1.0, 2.0));
    CHECK(x1(-1, -1) == std::complex<double>(1.0, 2.0));

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

    GLC x2 = GLC(2);
    x2.unserialize({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0});

    // In bounds
    // CHECK(x2(0) == 1.0);
    // CHECK(x2(1) == 2.0);
    // CHECK(x2(2) == 3.0);
    // CHECK(x2(3) == 4.0);
    // CHECK(x2(4) == 5.0);
    // CHECK(x2(5) == 6.0);
    // CHECK(x2(6) == 7.0);
    // CHECK(x2(7) == 8.0);

    // Out of bounds
    // CHECK(std::isnan(x2(-1)));
    // CHECK(std::isnan(x2(8)));

    // In bounds
    CHECK(x2(0, 0) == std::complex<double>(1.0, 2.0));
    CHECK(x2(0, 1) == std::complex<double>(3.0, 4.0));
    CHECK(x2(1, 0) == std::complex<double>(5.0, 6.0));
    CHECK(x2(1, 1) == std::complex<double>(7.0, 8.0));
    CHECK(x2(-1, -1) == std::complex<double>(7.0, 8.0));
    CHECK(x2(-1, -2) == std::complex<double>(5.0, 6.0));
    CHECK(x2(-2, -1) == std::complex<double>(3.0, 4.0));
    CHECK(x2(-2, -2) == std::complex<double>(1.0, 2.0));

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

TEST_CASE("GLC math_ops_GLC", "[domain]")
{
    using namespace Lielab::domain;

    GLC x1(2), x2(2);
    x1.unserialize({1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0});
    x2.unserialize({1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0});

    const GLC x1_prod_x2 = x1*x2;
    const Eigen::VectorXd x1_prod_x2bar = x1_prod_x2.serialize();
    REQUIRE(x1_prod_x2bar.size() == 8);
    CHECK(x1_prod_x2bar(0) == -4.75);
    CHECK(x1_prod_x2bar(1) == 12.0);
    CHECK(x1_prod_x2bar(2) == 0.75);
    CHECK(x1_prod_x2bar(3) == 15.5);
    CHECK(x1_prod_x2bar(4) == -1.25);
    CHECK(x1_prod_x2bar(5) == 5.75);
    CHECK(x1_prod_x2bar(6) == 2.75);
    CHECK(x1_prod_x2bar(7) == 6.75);

    x1 *= x2;
    const Eigen::VectorXd x1bar = x1.serialize();
    REQUIRE(x1bar.size() == 8);
    CHECK(x1bar(0) == -4.75);
    CHECK(x1bar(1) == 12.0);
    CHECK(x1bar(2) == 0.75);
    CHECK(x1bar(3) == 15.5);
    CHECK(x1bar(4) == -1.25);
    CHECK(x1bar(5) == 5.75);
    CHECK(x1bar(6) == 2.75);
    CHECK(x1bar(7) == 6.75);

    x1.unserialize({1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0});

    const GLC x1_inv = x1.inverse();
    const Eigen::VectorXd x1_invbar = x1_inv.serialize();
    REQUIRE(x1_invbar.size() == 8);
    CHECK_THAT(x1_invbar(0), Catch::Matchers::WithinULP(-0.25, 2));
    CHECK_THAT(x1_invbar(1), Catch::Matchers::WithinULP(0.25, 2));
    CHECK_THAT(x1_invbar(2), Catch::Matchers::WithinULP(1.0, 2));
    CHECK_THAT(x1_invbar(3), Catch::Matchers::WithinULP(-0.75, 4));
    CHECK_THAT(x1_invbar(4), Catch::Matchers::WithinULP(0.25, 2));
    CHECK_THAT(x1_invbar(5), Catch::Matchers::WithinULP(-0.25, 2));
    CHECK_THAT(x1_invbar(6), Catch::Matchers::WithinULP(-0.5, 2));
    CHECK_THAT(x1_invbar(7), Catch::Matchers::WithinULP(0.25, 2));
}

TEST_CASE("GLC project", "[domain]")
{
    using namespace Lielab::domain;

    const Eigen::MatrixXcd rand_2_2 = Eigen::MatrixXcd::Random(2, 2);
    const Eigen::MatrixXcd proj_2_2 = GLC::project(rand_2_2);

    REQUIRE(proj_2_2.rows() == 2);
    REQUIRE(proj_2_2.cols() == 2);
    CHECK(proj_2_2(0, 0) == rand_2_2(0, 0));
    CHECK(proj_2_2(0, 1) == rand_2_2(0, 1));
    CHECK(proj_2_2(1, 0) == rand_2_2(1, 0));
    CHECK(proj_2_2(1, 1) == rand_2_2(1, 1));

    const Eigen::MatrixXcd rand_3_3 = Eigen::MatrixXcd::Random(3, 3);
    const Eigen::MatrixXcd proj_3_3 = GLC::project(rand_3_3);

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

    const Eigen::MatrixXcd rand_2_3 = Eigen::MatrixXcd::Random(2, 3);
    const Eigen::MatrixXcd proj_2_3 = GLC::project(rand_2_3);

    REQUIRE(proj_2_3.rows() == 2);
    REQUIRE(proj_2_3.cols() == 2);
    CHECK(proj_2_3(0, 0) == rand_2_3(0, 0));
    CHECK(proj_2_3(0, 1) == rand_2_3(0, 1));
    CHECK(proj_2_3(1, 0) == rand_2_3(1, 0));
    CHECK(proj_2_3(1, 1) == rand_2_3(1, 1));

    const Eigen::MatrixXcd rand_3_2 = Eigen::MatrixXcd::Random(3, 2);
    const Eigen::MatrixXcd proj_3_2 = GLC::project(rand_3_2);

    REQUIRE(proj_3_2.rows() == 2);
    REQUIRE(proj_3_2.cols() == 2);
    CHECK(proj_3_2(0, 0) == rand_3_2(0, 0));
    CHECK(proj_3_2(0, 1) == rand_3_2(0, 1));
    CHECK(proj_3_2(1, 0) == rand_3_2(1, 0));
    CHECK(proj_3_2(1, 1) == rand_3_2(1, 1));
}

