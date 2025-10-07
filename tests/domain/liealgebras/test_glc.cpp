#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

#include <iostream>

TEST_CASE("glc to_string", "[domain]")
{
    using namespace Lielab::domain;

    const glc x0 = glc(0);
    CHECK(x0.to_string() == "gl(0, C)");
    const glc x1 = glc(1);
    CHECK(x1.to_string() == "gl(1, C)");
    const glc x10 = glc(10);
    CHECK(x10.to_string() == "gl(10, C)");
}

TEST_CASE("glc main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const glc xblank = glc();
    CHECK(xblank.get_dimension() == 0);

    const glc x0 = glc(0);
    CHECK(x0.get_dimension() == 0);
    const glc x1 = glc(1);
    CHECK(x1.get_dimension() == 2);
    const glc x10 = glc(10);
    CHECK(x10.get_dimension() == 200);
}

TEST_CASE("glc matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const glc x0 = glc(Eigen::MatrixXcd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    const glc x1 = glc(Eigen::MatrixXcd::Random(1, 1));
    CHECK(x1.get_shape() == 1);

    const glc x2 = glc(Eigen::MatrixXcd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(glc(Eigen::MatrixXcd::Random(2, 3)));
    CHECK_THROWS(glc(Eigen::MatrixXcd::Random(3, 2)));
}

TEST_CASE("glc basis_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const glc xm10 = glc::basis(-1, 0);
    CHECK(xm10.get_dimension() == 0);
    const Eigen::VectorXd xm10bar = xm10.get_vector();
    CHECK(xm10bar.size() == 0);

    const glc x00 = glc::basis(0, 0);
    CHECK(x00.get_dimension() == 0);
    const Eigen::VectorXd x00bar = x00.get_vector();
    CHECK(x00bar.size() == 0);

    const glc x10 = glc::basis(1, 0);
    CHECK(x10.get_dimension() == 0);
    const Eigen::VectorXd x10bar = x10.get_vector();
    CHECK(x10bar.size() == 0);

    const glc x01 = glc::basis(0, 1);
    CHECK(x01.get_dimension() == 2);
    const Eigen::VectorXd x01bar = x01.get_vector();
    REQUIRE(x01bar.size() == 2);
    CHECK(x01bar(0) == 1.0);
    CHECK(x01bar(1) == 0.0);

    const glc x11 = glc::basis(1, 1);
    CHECK(x11.get_dimension() == 2);
    const Eigen::VectorXd x11bar = x11.get_vector();
    REQUIRE(x11bar.size() == 2);
    CHECK(x11bar(0) == 0.0);
    CHECK(x11bar(1) == 1.0);

    const glc x21 = glc::basis(2, 1);
    CHECK(x21.get_dimension() == 2);
    const Eigen::VectorXd x21bar = x21.get_vector();
    REQUIRE(x21bar.size() == 2);
    CHECK(x21bar(0) == 0.0);
    CHECK(x21bar(1) == 0.0);

    const glc x02 = glc::basis(0, 2);
    CHECK(x02.get_dimension() == 8);
    const Eigen::VectorXd x02bar = x02.get_vector();
    REQUIRE(x02bar.size() == 8);
    CHECK(x02bar(0) == 1.0);
    CHECK(x02bar(1) == 0.0);
    CHECK(x02bar(2) == 0.0);
    CHECK(x02bar(3) == 0.0);
    CHECK(x02bar(4) == 0.0);
    CHECK(x02bar(5) == 0.0);
    CHECK(x02bar(6) == 0.0);
    CHECK(x02bar(7) == 0.0);
}

TEST_CASE("glc from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const glc x0 = glc::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXcd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const glc x1 = glc::from_shape(1);
    CHECK(x1.get_dimension() == 2);
    const Eigen::MatrixXcd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const glc x2 = glc::from_shape(2);
    CHECK(x2.get_dimension() == 8);
    const Eigen::MatrixXcd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("glc get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    glc zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
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

TEST_CASE("glc set/get_vector", "[domain]")
{
    /*!
    * Tests the set/get_vector operation.
    */

    using namespace Lielab::domain;

    glc x0 = glc(0);
    x0.set_vector({});
    Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0bar.size() == 0);

    glc x1 = glc(1);
    x1.set_vector({1.0, 2.0, 3.0});
    Eigen::VectorXd x1bar = x1.get_vector();

    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 1.0);
    CHECK(x1bar(1) == 2.0);

    x1.set_vector({4.0, 5.0});
    x1bar = x1.get_vector();

    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 4.0);
    CHECK(x1bar(1) == 5.0);

    x1.set_vector({6.0});
    x1bar = x1.get_vector();

    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 6.0);
    CHECK(x1bar(1) == 5.0);

    glc x2 = glc(2);
    x2.set_vector({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});
    Eigen::VectorXd x2bar = x2.get_vector();

    REQUIRE(x2bar.size() == 8);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 3.0);
    CHECK(x2bar(3) == 4.0);
    CHECK(x2bar(4) == 5.0);
    CHECK(x2bar(5) == 6.0);
    CHECK(x2bar(6) == 7.0);
    CHECK(x2bar(7) == 0.0);

    x2.set_vector({8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0});
    x2bar = x2.get_vector();

    REQUIRE(x2bar.size() == 8);
    CHECK(x2bar(0) == 8.0);
    CHECK(x2bar(1) == 9.0);
    CHECK(x2bar(2) == 10.0);
    CHECK(x2bar(3) == 11.0);
    CHECK(x2bar(4) == 12.0);
    CHECK(x2bar(5) == 13.0);
    CHECK(x2bar(6) == 14.0);
    CHECK(x2bar(7) == 15.0);

    x2.set_vector({16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0});
    x2bar = x2.get_vector();

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

TEST_CASE("glc get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    glc x0 = glc(0);
    x0.set_vector({});
    Eigen::MatrixXcd x0hat = x0.get_matrix();

    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    glc x1 = glc(1);
    x1.set_vector({1.0, 2.0, 3.0});
    Eigen::MatrixXcd x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == std::complex<double>(1.0, 2.0));

    x1.set_vector({4.0, 5.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == std::complex<double>(4.0, 5.0));

    x1.set_vector({6.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 1);
    REQUIRE(x1hat.cols() == 1);
    CHECK(x1hat(0, 0) == std::complex<double>(6.0, 5.0));

    glc x2 = glc(2);
    x2.set_vector({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});
    Eigen::MatrixXcd x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == std::complex<double>(1.0, 2.0));
    CHECK(x2hat(0, 1) == std::complex<double>(3.0, 4.0));
    CHECK(x2hat(1, 0) == std::complex<double>(5.0, 6.0));
    CHECK(x2hat(1, 1) == std::complex<double>(7.0, 0.0));

    x2.set_vector({8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == std::complex<double>(8.0, 9.0));
    CHECK(x2hat(0, 1) == std::complex<double>(10.0, 11.0));
    CHECK(x2hat(1, 0) == std::complex<double>(12.0, 13.0));
    CHECK(x2hat(1, 1) == std::complex<double>(14.0, 15.0));

    x2.set_vector({16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 2);
    REQUIRE(x2hat.cols() == 2);
    CHECK(x2hat(0, 0) == std::complex<double>(16.0, 17.0));
    CHECK(x2hat(0, 1) == std::complex<double>(18.0, 19.0));
    CHECK(x2hat(1, 0) == std::complex<double>(20.0, 21.0));
    CHECK(x2hat(1, 1) == std::complex<double>(22.0, 23.0));
}

TEST_CASE("glc operator()", "[domain]")
{
    using namespace Lielab::domain;

    glc x0 = glc::from_shape(0);
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

    glc x1 = glc(1);
    x1.set_vector({1.0, 2.0});

    // In bounds
    CHECK(x1(0) == 1.0);
    CHECK(x1(1) == 2.0);
    CHECK(x1(-1) == 2.0);
    CHECK(x1(-2) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x1(-3)));
    CHECK(std::isnan(x1(2)));

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

    glc x2 = glc(2);
    x2.set_vector({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0});

    // In bounds
    CHECK(x2(0) == 1.0);
    CHECK(x2(1) == 2.0);
    CHECK(x2(2) == 3.0);
    CHECK(x2(3) == 4.0);
    CHECK(x2(4) == 5.0);
    CHECK(x2(5) == 6.0);
    CHECK(x2(6) == 7.0);
    CHECK(x2(7) == 8.0);
    CHECK(x2(-1) == 8.0);
    CHECK(x2(-2) == 7.0);
    CHECK(x2(-3) == 6.0);
    CHECK(x2(-4) == 5.0);
    CHECK(x2(-5) == 4.0);
    CHECK(x2(-6) == 3.0);
    CHECK(x2(-7) == 2.0);
    CHECK(x2(-8) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x2(-9)));
    CHECK(std::isnan(x2(8)));

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

// TODO: mathops int here

TEST_CASE("glc math_ops_double", "[domain]")
{
    using namespace Lielab::domain;
    
    const std::complex<double> j(0.0, 1.0);

    glc x1(2);
    x1.set_vector({1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0});

    const glc x1_lm_2 = 2.0*x1;
    CHECK(x1_lm_2(0, 0) == std::complex<double>(2.5, 5.0));
    CHECK(x1_lm_2(0, 1) == std::complex<double>(7.5, 2.0));
    CHECK(x1_lm_2(1, 0) == std::complex<double>(2.0, 2.0));
    CHECK(x1_lm_2(1, 1) == std::complex<double>(2.0, 2.0));

    const glc x1_rm_2 = x1*2.0;
    CHECK(x1_rm_2(0, 0) == std::complex<double>(2.5, 5.0));
    CHECK(x1_rm_2(0, 1) == std::complex<double>(7.5, 2.0));
    CHECK(x1_rm_2(1, 0) == std::complex<double>(2.0, 2.0));
    CHECK(x1_rm_2(1, 1) == std::complex<double>(2.0, 2.0));

    x1 *= 2.0;
    CHECK(x1(0, 0) == std::complex<double>(2.5, 5.0));
    CHECK(x1(0, 1) == std::complex<double>(7.5, 2.0));
    CHECK(x1(1, 0) == std::complex<double>(2.0, 2.0));
    CHECK(x1(1, 1) == std::complex<double>(2.0, 2.0));

    x1.set_vector({1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0});

    const glc x1_lm_2j = (2.0*j)*x1;
    CHECK(x1_lm_2j(0, 0) == std::complex<double>(-5.0, 2.5));
    CHECK(x1_lm_2j(0, 1) == std::complex<double>(-2.0, 7.5));
    CHECK(x1_lm_2j(1, 0) == std::complex<double>(-2.0, 2.0));
    CHECK(x1_lm_2j(1, 1) == std::complex<double>(-2.0, 2.0));

    const glc x1_rm_2j = x1*(2.0*j);
    CHECK(x1_rm_2j(0, 0) == std::complex<double>(-5.0, 2.5));
    CHECK(x1_rm_2j(0, 1) == std::complex<double>(-2.0, 7.5));
    CHECK(x1_rm_2j(1, 0) == std::complex<double>(-2.0, 2.0));
    CHECK(x1_rm_2j(1, 1) == std::complex<double>(-2.0, 2.0));

    x1 *= 2.0*j;
    CHECK(x1(0, 0) == std::complex<double>(-5.0, 2.5));
    CHECK(x1(0, 1) == std::complex<double>(-2.0, 7.5));
    CHECK(x1(1, 0) == std::complex<double>(-2.0, 2.0));
    CHECK(x1(1, 1) == std::complex<double>(-2.0, 2.0));

    x1.set_vector({1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0});

    const glc x1_d_2 = x1/2.0;
    CHECK(x1_d_2(0, 0) == std::complex<double>(0.625, 1.25));
    CHECK(x1_d_2(0, 1) == std::complex<double>(1.875, 0.5));
    CHECK(x1_d_2(1, 0) == std::complex<double>(0.5, 0.5));
    CHECK(x1_d_2(1, 1) == std::complex<double>(0.5, 0.5));

    x1 /= 2.0;
    CHECK(x1(0, 0) == std::complex<double>(0.625, 1.25));
    CHECK(x1(0, 1) == std::complex<double>(1.875, 0.5));
    CHECK(x1(1, 0) == std::complex<double>(0.5, 0.5));
    CHECK(x1(1, 1) == std::complex<double>(0.5, 0.5));

    x1.set_vector({1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0});

    const glc x1_d_2j = x1/(2.0*j);
    CHECK(x1_d_2j(0, 0) == std::complex<double>(1.25, -0.625));
    CHECK(x1_d_2j(0, 1) == std::complex<double>(0.5, -1.875));
    CHECK(x1_d_2j(1, 0) == std::complex<double>(0.5, -0.5));
    CHECK(x1_d_2j(1, 1) == std::complex<double>(0.5, -0.5));

    x1 /= 2.0*j;
    CHECK(x1(0, 0) == std::complex<double>(1.25, -0.625));
    CHECK(x1(0, 1) == std::complex<double>(0.5, -1.875));
    CHECK(x1(1, 0) == std::complex<double>(0.5, -0.5));
    CHECK(x1(1, 1) == std::complex<double>(0.5, -0.5));
}

TEST_CASE("glc math_ops_cn", "[domain]")
{
    using namespace Lielab::domain;

    glc x1(2), x2(2);
    x1.set_vector({1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0});
    x2.set_vector({1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0});

    const glc x1_add_x2 = x1 + x2;
    CHECK(x1_add_x2(0, 0) == std::complex<double>(2.25, 4.5));
    CHECK(x1_add_x2(0, 1) == std::complex<double>(6.75, 5.0));
    CHECK(x1_add_x2(1, 0) == std::complex<double>(2.0, 2.0));
    CHECK(x1_add_x2(1, 1) == std::complex<double>(2.0, 2.0));

    x1 += x2;
    CHECK(x1(0, 0) == std::complex<double>(2.25, 4.5));
    CHECK(x1(0, 1) == std::complex<double>(6.75, 5.0));
    CHECK(x1(1, 0) == std::complex<double>(2.0, 2.0));
    CHECK(x1(1, 1) == std::complex<double>(2.0, 2.0));

    x1.set_vector({1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0});

    const glc x1_sub_x2 = x1 - x2;
    CHECK(x1_sub_x2(0, 0) == std::complex<double>(-0.25, -0.5));
    CHECK(x1_sub_x2(0, 1) == std::complex<double>(-0.75, 3.0));
    CHECK(x1_sub_x2(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x1_sub_x2(1, 1) == std::complex<double>(0.0, 0.0));

    x1 -= x2;
    CHECK(x1(0, 0) == std::complex<double>(-0.25, -0.5));
    CHECK(x1(0, 1) == std::complex<double>(-0.75, 3.0));
    CHECK(x1(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x1(1, 1) == std::complex<double>(0.0, 0.0));

    x1.set_vector({1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0});

    const glc x1_unary_sub = (-x1);
    CHECK(x1_unary_sub(0, 0) == std::complex<double>(-1.0, -2.0));
    CHECK(x1_unary_sub(0, 1) == std::complex<double>(-3.0, -4.0));
    CHECK(x1_unary_sub(1, 0) == std::complex<double>(-1.0, -1.0));
    CHECK(x1_unary_sub(1, 1) == std::complex<double>(-1.0, -1.0));
}

TEST_CASE("glc from_vector", "[domain]")
{
    /*!
    * Tests the from_vector operation.
    */

    using namespace Lielab::domain;

    const glc x0 = glc::from_vector({});
    const Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0.get_shape() == 0);
    CHECK(x0bar.size() == 0);

    const glc x1 = glc::from_vector({1.0});
    const Eigen::VectorXd x1bar = x1.get_vector();

    CHECK(x1.get_shape() == 1);
    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 1.0);
    CHECK(x1bar(1) == 0.0);

    const glc x2 = glc::from_vector({1.0, 2.0});
    const Eigen::VectorXd x2bar = x2.get_vector();

    CHECK(x2.get_shape() == 1);
    REQUIRE(x2bar.size() == 2);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);

    const glc x3 = glc::from_vector({1.0, 2.0, 3.0});
    const Eigen::VectorXd x3bar = x3.get_vector();

    CHECK(x3.get_shape() == 2);
    REQUIRE(x3bar.size() == 8);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
    CHECK(x3bar(3) == 0.0);
    CHECK(x3bar(4) == 0.0);
    CHECK(x3bar(5) == 0.0);
    CHECK(x3bar(6) == 0.0);
    CHECK(x3bar(7) == 0.0);
}

TEST_CASE("glc from_complex_vector", "[domain]")
{
    /*!
    * Tests the from_complex_vector operation.
    */

    using namespace Lielab::domain;
    const std::complex<double> j(0.0, 1.0);

    const glc x0 = glc::from_complex_vector({});
    const Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0.get_shape() == 0);
    CHECK(x0bar.size() == 0);

    const glc x1 = glc::from_complex_vector({1.0});
    const Eigen::VectorXd x1bar = x1.get_vector();

    CHECK(x1.get_shape() == 1);
    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 1.0);
    CHECK(x1bar(1) == 0.0);

    const glc x2 = glc::from_complex_vector({1.0 + 2.0*j});
    const Eigen::VectorXd x2bar = x2.get_vector();

    CHECK(x2.get_shape() == 1);
    REQUIRE(x2bar.size() == 2);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);

    const glc x3 = glc::from_complex_vector({1.0 + 2.0*j, 3.0});
    const Eigen::VectorXd x3bar = x3.get_vector();

    CHECK(x3.get_shape() == 2);
    REQUIRE(x3bar.size() == 8);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
    CHECK(x3bar(3) == 0.0);
    CHECK(x3bar(4) == 0.0);
    CHECK(x3bar(5) == 0.0);
    CHECK(x3bar(6) == 0.0);
    CHECK(x3bar(7) == 0.0);
}

TEST_CASE("glc project", "[domain]")
{
    using namespace Lielab::domain;

    const Eigen::MatrixXcd rand_2_2 = Eigen::MatrixXcd::Random(2, 2);
    const Eigen::MatrixXcd proj_2_2 = glc::project(rand_2_2);

    REQUIRE(proj_2_2.rows() == 2);
    REQUIRE(proj_2_2.cols() == 2);
    CHECK(proj_2_2(0, 0) == rand_2_2(0, 0));
    CHECK(proj_2_2(0, 1) == rand_2_2(0, 1));
    CHECK(proj_2_2(1, 0) == rand_2_2(1, 0));
    CHECK(proj_2_2(1, 1) == rand_2_2(1, 1));

    const Eigen::MatrixXcd rand_3_3 = Eigen::MatrixXcd::Random(3, 3);
    const Eigen::MatrixXcd proj_3_3 = glc::project(rand_3_3);

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
    const Eigen::MatrixXcd proj_2_3 = glc::project(rand_2_3);

    REQUIRE(proj_2_3.rows() == 2);
    REQUIRE(proj_2_3.cols() == 2);
    CHECK(proj_2_3(0, 0) == rand_2_3(0, 0));
    CHECK(proj_2_3(0, 1) == rand_2_3(0, 1));
    CHECK(proj_2_3(1, 0) == rand_2_3(1, 0));
    CHECK(proj_2_3(1, 1) == rand_2_3(1, 1));

    const Eigen::MatrixXcd rand_3_2 = Eigen::MatrixXcd::Random(3, 2);
    const Eigen::MatrixXcd proj_3_2 = glc::project(rand_3_2);

    REQUIRE(proj_3_2.rows() == 2);
    REQUIRE(proj_3_2.cols() == 2);
    CHECK(proj_3_2(0, 0) == rand_3_2(0, 0));
    CHECK(proj_3_2(0, 1) == rand_3_2(0, 1));
    CHECK(proj_3_2(1, 0) == rand_3_2(1, 0));
    CHECK(proj_3_2(1, 1) == rand_3_2(1, 1));
}

TEST_CASE("glc get/from_vector", "[domain]")
{
    /*!
    * Tests the get/from vector operation.
    */

    using namespace Lielab::domain;

    const glc x0 = glc::from_vector({});
    const Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0bar.size() == 0);
    CHECK(x0.get_dimension() == 0);
    
    const glc x1 = glc::from_vector({1.0});
    const Eigen::VectorXd x1bar = x1.get_vector();

    CHECK(x1bar.size() == 2);
    CHECK(x1.get_dimension() == 2);
    CHECK(x1bar(0) == 1.0);
    CHECK(x1bar(1) == 0.0);

    const glc x2 = glc::from_vector({1.0, 2.0});
    const Eigen::VectorXd x2bar = x2.get_vector();

    CHECK(x2bar.size() == 2);
    CHECK(x2.get_dimension() == 2);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);

    const glc x3 = glc::from_vector({1.0, 2.0, 3.0});
    const Eigen::VectorXd x3bar = x3.get_vector();

    CHECK(x3bar.size() == 8);
    CHECK(x3.get_dimension() == 8);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
    CHECK(x3bar(3) == 0.0);
    CHECK(x3bar(4) == 0.0);
    CHECK(x3bar(5) == 0.0);
    CHECK(x3bar(6) == 0.0);
    CHECK(x3bar(7) == 0.0);

    const glc x4 = glc::from_vector({1.0, 2.0, 3.0, 4.0});
    const Eigen::VectorXd x4bar = x4.get_vector();

    CHECK(x4bar.size() == 8);
    CHECK(x4.get_dimension() == 8);
    CHECK(x4bar(0) == 1.0);
    CHECK(x4bar(1) == 2.0);
    CHECK(x4bar(2) == 3.0);
    CHECK(x4bar(3) == 4.0);
    CHECK(x4bar(4) == 0.0);
    CHECK(x4bar(5) == 0.0);
    CHECK(x4bar(6) == 0.0);
    CHECK(x4bar(7) == 0.0);

    const glc x5 = glc::from_vector({1.0, 2.0, 3.0, 4.0, 5.0});
    const Eigen::VectorXd x5bar = x5.get_vector();

    CHECK(x5bar.size() == 8);
    CHECK(x5.get_dimension() == 8);
    CHECK(x5bar(0) == 1.0);
    CHECK(x5bar(1) == 2.0);
    CHECK(x5bar(2) == 3.0);
    CHECK(x5bar(3) == 4.0);
    CHECK(x5bar(4) == 5.0);
    CHECK(x5bar(5) == 0.0);
    CHECK(x5bar(6) == 0.0);
    CHECK(x5bar(7) == 0.0);
}
