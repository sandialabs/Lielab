#include <Lielab.hpp>
#include <iostream>

#include <catch2/catch_all.hpp>

TEST_CASE("CN to_string", "[domain]")
{
    using namespace Lielab::domain;

    const CN xzero = CN::from_shape(0);
    CHECK(xzero.to_string() == "C^nan");
    const CN x0 = CN(0);
    CHECK(x0.to_string() == "C^0");
    const CN x1 = CN(1);
    CHECK(x1.to_string() == "C^1");
    const CN x10 = CN(10);
    CHECK(x10.to_string() == "C^10");
}

TEST_CASE("CN main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CN xblank = CN();
    CHECK(xblank.get_dimension() == 0);

    const CN x0 = CN(0);
    CHECK(x0.get_dimension() == 0);
    const CN x1 = CN(1);
    CHECK(x1.get_dimension() == 2);
    const CN x10 = CN(10);
    CHECK(x10.get_dimension() == 20);
}

TEST_CASE("CN matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CN x0 = CN(Eigen::MatrixXcd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    const CN x1 = CN(Eigen::MatrixXcd::Random(1, 1));
    CHECK(x1.get_shape() == 1);

    const CN x2 = CN(Eigen::MatrixXcd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(CN(Eigen::MatrixXcd::Random(2, 3)));
    CHECK_THROWS(CN(Eigen::MatrixXcd::Random(3, 2)));
}

TEST_CASE("CN from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CN x0 = CN::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXcd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const CN x1 = CN::from_shape(1);
    CHECK(x1.get_dimension() == 0);
    const Eigen::MatrixXcd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const CN x2 = CN::from_shape(2);
    CHECK(x2.get_dimension() == 2);
    const Eigen::MatrixXcd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("CN get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    CN veryzero = CN::from_shape(0);
    CN zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    CHECK(veryzero.get_dimension() == 0);
    CHECK(zero.get_dimension() == 0);
    CHECK(one.get_dimension() == 2);
    CHECK(two.get_dimension() == 4);
    CHECK(three.get_dimension() == 6);
    CHECK(four.get_dimension() == 8);
    CHECK(five.get_dimension() == 10);
    CHECK(six.get_dimension() == 12);
    CHECK(seven.get_dimension() == 14);
    CHECK(eight.get_dimension() == 16);
}

TEST_CASE("CN serialize/unserialize", "[domain]")
{
    /*!
    * Tests the serialize/unserialize operation.
    */

    using namespace Lielab::domain;

    CN xzero = CN::from_shape(0);
    xzero.unserialize({});
    Eigen::VectorXd xzerobar = xzero.serialize();

    CHECK(xzerobar.size() == 0);

    CN x0 = CN(0);
    x0.unserialize({});
    Eigen::VectorXd x0bar = x0.serialize();

    CHECK(x0bar.size() == 0);

    CN x1 = CN(1);
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

    CN x2 = CN(2);
    x2.unserialize({1.0, 2.0, 3.0});
    Eigen::VectorXd x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 3.0);
    CHECK(x2bar(3) == 0.0);

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

TEST_CASE("CN get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    CN xzero = CN::from_shape(0);
    xzero.unserialize({});
    Eigen::MatrixXcd xzerohat = xzero.get_matrix();

    CHECK(xzerohat.rows() == 0);
    CHECK(xzerohat.cols() == 0);

    CN x0 = CN(0);
    x0.unserialize({});
    Eigen::MatrixXcd x0hat = x0.get_matrix();

    REQUIRE(x0hat.rows() == 1);
    REQUIRE(x0hat.cols() == 1);
    CHECK(x0hat(0, 0) == std::complex<double>(1.0, 0.0));

    CN x1 = CN(1);
    x1.unserialize({1.0, 2.0, 3.0});
    Eigen::MatrixXcd x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 1) == std::complex<double>(1.0, 2.0));

    CHECK(x1hat(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(x1hat(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x1hat(1, 1) == std::complex<double>(1.0, 0.0));

    x1.unserialize({4.0, 5.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 1) == std::complex<double>(4.0, 5.0));

    CHECK(x1hat(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(x1hat(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x1hat(1, 1) == std::complex<double>(1.0, 0.0));

    x1.unserialize({6.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 1) == std::complex<double>(6.0, 5.0));

    CHECK(x1hat(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(x1hat(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x1hat(1, 1) == std::complex<double>(1.0, 0.0));

    CN x2 = CN(2);
    x2.unserialize({1.0, 2.0, 3.0});
    Eigen::MatrixXcd x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 2) == std::complex<double>(1.0, 2.0));
    CHECK(x2hat(1, 2) == std::complex<double>(3.0, 0.0));

    CHECK(x2hat(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(x2hat(0, 1) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(1, 1) == std::complex<double>(1.0, 0.0));
    CHECK(x2hat(2, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(2, 1) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(2, 2) == std::complex<double>(1.0, 0.0));

    x2.unserialize({4.0, 5.0, 6.0, 7.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 2) == std::complex<double>(4.0, 5.0));
    CHECK(x2hat(1, 2) == std::complex<double>(6.0, 7.0));

    CHECK(x2hat(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(x2hat(0, 1) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(1, 1) == std::complex<double>(1.0, 0.0));
    CHECK(x2hat(2, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(2, 1) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(2, 2) == std::complex<double>(1.0, 0.0));

    x2.unserialize({8.0, 9.0, 10.0, 11.0, 12.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 2) == std::complex<double>(8.0, 9.0));
    CHECK(x2hat(1, 2) == std::complex<double>(10.0, 11.0));

    CHECK(x2hat(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(x2hat(0, 1) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(1, 1) == std::complex<double>(1.0, 0.0));
    CHECK(x2hat(2, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(2, 1) == std::complex<double>(0.0, 0.0));
    CHECK(x2hat(2, 2) == std::complex<double>(1.0, 0.0));
}

TEST_CASE("CN operator()", "[domain]")
{
    using namespace Lielab::domain;

    CN xzero = CN::from_shape(0);
    xzero.unserialize({});

    // Out of bounds
    CHECK(std::isnan(xzero(-1)));
    CHECK(std::isnan(xzero(0)));
    CHECK(std::isnan(xzero(1)));

    // Out of bounds
    CHECK(std::isnan(xzero(0, -1).real()));
    CHECK(std::isnan(xzero(0, -1).imag()));
    CHECK(std::isnan(xzero(-1, 0).real()));
    CHECK(std::isnan(xzero(-1, 0).imag()));
    CHECK(std::isnan(xzero(-1, -1).real()));
    CHECK(std::isnan(xzero(-1, -1).imag()));
    CHECK(std::isnan(xzero(0, 0).real()));
    CHECK(std::isnan(xzero(0, 0).imag()));
    CHECK(std::isnan(xzero(0, 1).real()));
    CHECK(std::isnan(xzero(0, 1).imag()));
    CHECK(std::isnan(xzero(1, 0).real()));
    CHECK(std::isnan(xzero(1, 0).imag()));
    CHECK(std::isnan(xzero(1, 1).real()));
    CHECK(std::isnan(xzero(1, 1).imag()));

    CN x1 = CN(1);
    x1.unserialize({1.0, 2.0});

    // In bounds
    CHECK(x1(0) == 1.0);
    CHECK(x1(1) == 2.0);
    CHECK(x1(-1) == 2.0);
    CHECK(x1(-2) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x1(-3)));
    CHECK(std::isnan(x1(2)));

    // In bounds
    CHECK(x1(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(x1(0, 1) == std::complex<double>(1.0, 2.0));
    CHECK(x1(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x1(1, 1) == std::complex<double>(1.0, 0.0));
    CHECK(x1(-1, -1) == std::complex<double>(1.0, 0.0));
    CHECK(x1(-1, -2) == std::complex<double>(0.0, 0.0));
    CHECK(x1(-2, -1) == std::complex<double>(1.0, 2.0));
    CHECK(x1(-2, -2) == std::complex<double>(1.0, 0.0));

    // Out of bounds
    CHECK(std::isnan(x1(0, -3).real()));
    CHECK(std::isnan(x1(0, -3).imag()));
    CHECK(std::isnan(x1(-3, 0).real()));
    CHECK(std::isnan(x1(-3, 0).imag()));
    CHECK(std::isnan(x1(-3, -3).real()));
    CHECK(std::isnan(x1(-3, -3).imag()));
    CHECK(std::isnan(x1(0, 2).real()));
    CHECK(std::isnan(x1(0, 2).imag()));
    CHECK(std::isnan(x1(2, 0).real()));
    CHECK(std::isnan(x1(2, 0).imag()));
    CHECK(std::isnan(x1(2, 2).real()));
    CHECK(std::isnan(x1(2, 2).imag()));

    CN x2 = CN(2);
    x2.unserialize({1.0, 2.0, 3.0, 4.0});

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
    CHECK(x2(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(x2(0, 1) == std::complex<double>(0.0, 0.0));
    CHECK(x2(0, 2) == std::complex<double>(1.0, 2.0));
    CHECK(x2(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x2(1, 1) == std::complex<double>(1.0, 0.0));
    CHECK(x2(1, 2) == std::complex<double>(3.0, 4.0));
    CHECK(x2(2, 0) == std::complex<double>(0.0, 0.0));
    CHECK(x2(2, 1) == std::complex<double>(0.0, 0.0));
    CHECK(x2(2, 2) == std::complex<double>(1.0, 0.0));
    CHECK(x2(-1, -1) == std::complex<double>(1.0, 0.0));
    CHECK(x2(-1, -2) == std::complex<double>(0.0, 0.0));
    CHECK(x2(-1, -3) == std::complex<double>(0.0, 0.0));
    CHECK(x2(-2, -1) == std::complex<double>(3.0, 4.0));
    CHECK(x2(-2, -2) == std::complex<double>(1.0, 0.0));
    CHECK(x2(-2, -3) == std::complex<double>(0.0, 0.0));
    CHECK(x2(-3, -1) == std::complex<double>(1.0, 2.0));
    CHECK(x2(-3, -2) == std::complex<double>(0.0, 0.0));
    CHECK(x2(-3, -3) == std::complex<double>(1.0, 0.0));

    // Out of bounds
    CHECK(std::isnan(x2(0, -4).real()));
    CHECK(std::isnan(x2(0, -4).imag()));
    CHECK(std::isnan(x2(-4, 0).real()));
    CHECK(std::isnan(x2(-4, 0).imag()));
    CHECK(std::isnan(x2(-4, -4).real()));
    CHECK(std::isnan(x2(-4, -4).imag()));
    CHECK(std::isnan(x2(0, 3).real()));
    CHECK(std::isnan(x2(0, 3).imag()));
    CHECK(std::isnan(x2(3, 0).real()));
    CHECK(std::isnan(x2(3, 0).imag()));
    CHECK(std::isnan(x2(3, 3).real()));
    CHECK(std::isnan(x2(3, 3).imag()));
}

TEST_CASE("CN operator[]", "[domain]")
{
    using namespace Lielab::domain;

    CN xzero = CN::from_shape(0);
    xzero.unserialize({});

    // Out of bounds
    CHECK(std::isnan(xzero[-1].real()));
    CHECK(std::isnan(xzero[-1].imag()));
    CHECK(std::isnan(xzero[0].real()));
    CHECK(std::isnan(xzero[0].imag()));
    CHECK(std::isnan(xzero[1].real()));
    CHECK(std::isnan(xzero[1].imag()));

    CN x1 = CN(1);
    x1.unserialize({1.0, 2.0});

    // In bounds
    CHECK(x1[0] == std::complex<double>(1.0, 2.0));
    CHECK(x1[-1] == std::complex<double>(1.0, 2.0));

    // Out of bounds
    CHECK(std::isnan(x1[-2].real()));
    CHECK(std::isnan(x1[-2].imag()));
    CHECK(std::isnan(x1[1].real()));
    CHECK(std::isnan(x1[1].imag()));

    CN x2 = CN(2);
    x2.unserialize({1.0, 2.0, 3.0, 4.0});

    // In bounds
    CHECK(x2[0] == std::complex<double>(1.0, 2.0));
    CHECK(x2[1] == std::complex<double>(3.0, 4.0));
    CHECK(x2[-1] == std::complex<double>(3.0, 4.0));
    CHECK(x2[-2] == std::complex<double>(1.0, 2.0));

    // Out of bounds
    CHECK(std::isnan(x2[-3].real()));
    CHECK(std::isnan(x2[-3].imag()));
    CHECK(std::isnan(x2[2].real()));
    CHECK(std::isnan(x2[2].imag()));
}

TEST_CASE("CN math_ops_CN", "[domain]")
{
    using namespace Lielab::domain;

    CN x1(2), x2(2);
    x1.unserialize({1.0, 2.0, 3.0, 4.0});
    x2.unserialize({1.25, 2.5, 3.75, 1.0});

    const CN x1_prod_x2 = x1*x2;
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

    const CN x1_inv = x1.inverse();
    const Eigen::VectorXd x1_invbar = x1_inv.serialize();
    REQUIRE(x1_invbar.size() == 4);
    CHECK(x1_invbar(0) == -1.0);
    CHECK(x1_invbar(1) == -2.0);
    CHECK(x1_invbar(2) == -3.0);
    CHECK(x1_invbar(3) == -4.0);
}

TEST_CASE("CN from_vector", "[domain]")
{
    /*!
    * Tests the from_vector operation.
    */

    using namespace Lielab::domain;

    const CN x0 = CN::from_vector({});
    const Eigen::VectorXd x0bar = x0.serialize();

    CHECK(x0.get_shape() == 1);
    CHECK(x0bar.size() == 0);

    const CN x1 = CN::from_vector({1.0});
    const Eigen::VectorXd x1bar = x1.serialize();

    CHECK(x1.get_shape() == 2);
    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 1.0);
    CHECK(x1bar(1) == 0.0);

    const CN x2 = CN::from_vector({1.0, 2.0});
    const Eigen::VectorXd x2bar = x2.serialize();

    CHECK(x2.get_shape() == 2);
    REQUIRE(x2bar.size() == 2);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);

    const CN x3 = CN::from_vector({1.0, 2.0, 3.0});
    const Eigen::VectorXd x3bar = x3.serialize();

    CHECK(x3.get_shape() == 3);
    REQUIRE(x3bar.size() == 4);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
    CHECK(x3bar(3) == 0.0);
}

TEST_CASE("CN from/to_complex_vector", "[domain]")
{
    /*!
    * Tests the from/to_complex_vector operations.
    */

    using namespace Lielab::domain;
    const std::complex<double> j(0.0, 1.0);

    const CN x0 = CN::from_complex_vector({});
    const Eigen::VectorXcd x0bar = x0.to_complex_vector();

    CHECK(x0.get_shape() == 1);
    CHECK(x0bar.size() == 0);

    const CN x1 = CN::from_complex_vector({1.0});
    const Eigen::VectorXcd x1bar = x1.to_complex_vector();

    CHECK(x1.get_shape() == 2);
    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == std::complex<double>(1.0, 0.0));

    const CN x2 = CN::from_complex_vector({1.0 + 2.0*j});
    const Eigen::VectorXcd x2bar = x2.to_complex_vector();

    CHECK(x2.get_shape() == 2);
    REQUIRE(x2bar.size() == 1);
    CHECK(x2bar(0) == std::complex<double>(1.0, 2.0));

    const CN x3 = CN::from_complex_vector({1.0 + 2.0*j, 3.0});
    const Eigen::VectorXcd x3bar = x3.to_complex_vector();

    CHECK(x3.get_shape() == 3);
    REQUIRE(x3bar.size() == 2);
    CHECK(x3bar(0) == std::complex<double>(1.0, 2.0));
    CHECK(x3bar(1) == std::complex<double>(3.0, 0.0));
}

TEST_CASE("CN project", "[domain]")
{
    using namespace Lielab::domain;

    const Eigen::MatrixXcd rand_2_2 = Eigen::MatrixXcd::Random(2, 2);
    const Eigen::MatrixXcd proj_2_2 = CN::project(rand_2_2);

    REQUIRE(proj_2_2.rows() == 2);
    REQUIRE(proj_2_2.cols() == 2);
    CHECK(proj_2_2(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(proj_2_2(0, 1) == rand_2_2(0, 1));
    CHECK(proj_2_2(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(proj_2_2(1, 1) == std::complex<double>(1.0, 0.0));

    const Eigen::MatrixXcd rand_3_3 = Eigen::MatrixXcd::Random(3, 3);
    const Eigen::MatrixXcd proj_3_3 = CN::project(rand_3_3);

    REQUIRE(proj_3_3.rows() == 3);
    REQUIRE(proj_3_3.cols() == 3);
    CHECK(proj_3_3(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(proj_3_3(0, 1) == std::complex<double>(0.0, 0.0));
    CHECK(proj_3_3(0, 2) == rand_3_3(0, 2));
    CHECK(proj_3_3(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(proj_3_3(1, 1) == std::complex<double>(1.0, 0.0));
    CHECK(proj_3_3(1, 2) == rand_3_3(1, 2));
    CHECK(proj_3_3(2, 0) == std::complex<double>(0.0, 0.0));
    CHECK(proj_3_3(2, 1) == std::complex<double>(0.0, 0.0));
    CHECK(proj_3_3(2, 2) == std::complex<double>(1.0, 0.0));

    const Eigen::MatrixXcd rand_2_3 = Eigen::MatrixXcd::Random(2, 3);
    const Eigen::MatrixXcd proj_2_3 = CN::project(rand_2_3);

    REQUIRE(proj_2_3.rows() == 2);
    REQUIRE(proj_2_3.cols() == 2);
    CHECK(proj_2_3(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(proj_2_3(0, 1) == rand_2_3(0, 1));
    CHECK(proj_2_3(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(proj_2_3(1, 1) == std::complex<double>(1.0, 0.0));

    const Eigen::MatrixXcd rand_3_2 = Eigen::MatrixXcd::Random(3, 2);
    const Eigen::MatrixXcd proj_3_2 = CN::project(rand_3_2);

    REQUIRE(proj_3_2.rows() == 2);
    REQUIRE(proj_3_2.cols() == 2);
    CHECK(proj_3_2(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(proj_3_2(0, 1) == rand_3_2(0, 1));
    CHECK(proj_3_2(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(proj_3_2(1, 1) == std::complex<double>(1.0, 0.0));
}
