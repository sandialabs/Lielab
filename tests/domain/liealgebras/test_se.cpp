#include <Lielab.hpp>
#include <iostream>

#include <catch2/catch_all.hpp>

TEST_CASE("se to_string", "[domain]")
{
    using namespace Lielab::domain;

    const se xzero = se::from_shape(0);
    CHECK(xzero.to_string() == "se(nan)");
    const se x0 = se(0);
    CHECK(x0.to_string() == "se(0)");
    const se x1 = se(1);
    CHECK(x1.to_string() == "se(1)");
    const se x10 = se(10);
    CHECK(x10.to_string() == "se(10)");
}

TEST_CASE("se main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const se xblank = se();
    CHECK(xblank.get_dimension() == 0);

    const se x0 = se(0);
    CHECK(x0.get_dimension() == 0);
    const se x1 = se(1);
    CHECK(x1.get_dimension() == 1);
    const se x3 = se(3);
    CHECK(x3.get_dimension() == 6);
}

TEST_CASE("se matrix_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const se x0 = se(Eigen::MatrixXd::Random(0, 0));
    CHECK(x0.get_shape() == 0);

    const se x1 = se(Eigen::MatrixXd::Random(1, 1));
    CHECK(x1.get_shape() == 1);

    const se x2 = se(Eigen::MatrixXd::Random(2, 2));
    CHECK(x2.get_shape() == 2);

    CHECK_THROWS(se(Eigen::MatrixXd::Random(2, 3)));
    CHECK_THROWS(se(Eigen::MatrixXd::Random(3, 2)));
}

TEST_CASE("se basis_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const se xm10 = se::basis(-1, 0);
    CHECK(xm10.get_dimension() == 0);
    const Eigen::VectorXd xm10bar = xm10.get_vector();
    CHECK(xm10bar.size() == 0);

    const se x00 = se::basis(0, 0);
    CHECK(x00.get_dimension() == 0);
    const Eigen::VectorXd x00bar = x00.get_vector();
    CHECK(x00bar.size() == 0);

    const se x10 = se::basis(1, 0);
    CHECK(x10.get_dimension() == 0);
    const Eigen::VectorXd x10bar = x10.get_vector();
    CHECK(x10bar.size() == 0);

    const se x01 = se::basis(0, 1);
    CHECK(x01.get_dimension() == 1);
    const Eigen::VectorXd x01bar = x01.get_vector();
    REQUIRE(x01bar.size() == 1);
    CHECK(x01bar(0) == 1.0);

    const se x11 = se::basis(1, 1);
    CHECK(x11.get_dimension() == 1);
    const Eigen::VectorXd x11bar = x11.get_vector();
    REQUIRE(x11bar.size() == 1);
    CHECK(x11bar(0) == 0.0);

    const se x21 = se::basis(2, 1);
    CHECK(x21.get_dimension() == 1);
    const Eigen::VectorXd x21bar = x21.get_vector();
    REQUIRE(x21bar.size() == 1);
    CHECK(x21bar(0) == 0.0);

    const se x02 = se::basis(0, 2);
    CHECK(x02.get_dimension() == 3);
    const Eigen::VectorXd x02bar = x02.get_vector();
    REQUIRE(x02bar.size() == 3);
    CHECK(x02bar(0) == 1.0);
    CHECK(x02bar(1) == 0.0);
    CHECK(x02bar(2) == 0.0);
}

TEST_CASE("se from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const se x0 = se::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const se x1 = se::from_shape(1);
    CHECK(x1.get_dimension() == 0);
    const Eigen::MatrixXd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const se x2 = se::from_shape(2);
    CHECK(x2.get_dimension() == 1);
    const Eigen::MatrixXd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("se get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    se veryzero = se::from_shape(0);
    se zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

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

TEST_CASE("se set/get_vector", "[domain]")
{
    /*!
    * Tests the set/get_vector operation.
    */

    using namespace Lielab::domain;

    se xzero = se::from_shape(0);
    xzero.set_vector({});
    Eigen::VectorXd xzerobar = xzero.get_vector();

    CHECK(xzerobar.size() == 0);

    se x0 = se(0);
    x0.set_vector({});
    Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0bar.size() == 0);

    se x2 = se(2);
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

    se x3 = se(3);
    x3.set_vector({1.0, 2.0, 3.0, 4.0, 5.0});
    Eigen::VectorXd x3bar = x3.get_vector();

    REQUIRE(x3bar.size() == 6);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
    CHECK(x3bar(3) == 4.0);
    CHECK(x3bar(4) == 5.0);
    CHECK(x3bar(5) == 0.0);

    x3.set_vector({7.0, 8.0, 9.0, 10.0, 11.0, 12.0});
    x3bar = x3.get_vector();

    REQUIRE(x3bar.size() == 6);
    CHECK(x3bar(0) == 7.0);
    CHECK(x3bar(1) == 8.0);
    CHECK(x3bar(2) == 9.0);
    CHECK(x3bar(3) == 10.0);
    CHECK(x3bar(4) == 11.0);
    CHECK(x3bar(5) == 12.0);

    x3.set_vector({13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0});
    x3bar = x3.get_vector();

    REQUIRE(x3bar.size() == 6);
    CHECK(x3bar(0) == 13.0);
    CHECK(x3bar(1) == 14.0);
    CHECK(x3bar(2) == 15.0);
    CHECK(x3bar(3) == 16.0);
    CHECK(x3bar(4) == 17.0);
    CHECK(x3bar(5) == 18.0);
}

TEST_CASE("se get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    se xzero = se::from_shape(0);
    xzero.set_vector({});
    Eigen::MatrixXd xzerohat = xzero.get_matrix();

    CHECK(xzerohat.rows() == 0);
    CHECK(xzerohat.cols() == 0);

    se x0 = se(0);
    x0.set_vector({});
    Eigen::MatrixXd x0hat = x0.get_matrix();

    REQUIRE(x0hat.rows() == 1);
    REQUIRE(x0hat.cols() == 1);
    CHECK(x0hat(0, 0) == 0.0);

    se x1 = se(1);
    x1.set_vector({1.0, 2.0});
    Eigen::MatrixXd x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 1) == 1.0);

    CHECK(x1hat(0, 0) == 0.0);
    CHECK(x1hat(1, 0) == 0.0);
    CHECK(x1hat(1, 1) == 0.0);

    x1.set_vector({3.0});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 1) == 3.0);

    CHECK(x1hat(0, 0) == 0.0);
    CHECK(x1hat(1, 0) == 0.0);
    CHECK(x1hat(1, 1) == 0.0);

    x1.set_vector({});
    x1hat = x1.get_matrix();

    REQUIRE(x1hat.rows() == 2);
    REQUIRE(x1hat.cols() == 2);
    CHECK(x1hat(0, 1) == 3.0);

    CHECK(x1hat(0, 0) == 0.0);
    CHECK(x1hat(1, 0) == 0.0);
    CHECK(x1hat(1, 1) == 0.0);

    se x2 = se(2);
    x2.set_vector({1.0});
    Eigen::MatrixXd x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 2) == 1.0);
    CHECK(x2hat(1, 2) == 0.0);

    CHECK(x2hat(0, 0) == 0.0);
    CHECK(x2hat(0, 1) == 0.0);
    CHECK(x2hat(1, 0) == 0.0);
    CHECK(x2hat(1, 1) == 0.0);
    CHECK(x2hat(2, 0) == 0.0);
    CHECK(x2hat(2, 1) == 0.0);
    CHECK(x2hat(2, 2) == 0.0);

    x2.set_vector({2.0, 3.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 2) == 2.0);
    CHECK(x2hat(1, 2) == 3.0);

    CHECK(x2hat(0, 0) == 0.0);
    CHECK(x2hat(0, 1) == 0.0);
    CHECK(x2hat(1, 0) == 0.0);
    CHECK(x2hat(1, 1) == 0.0);
    CHECK(x2hat(2, 0) == 0.0);
    CHECK(x2hat(2, 1) == 0.0);
    CHECK(x2hat(2, 2) == 0.0);

    const Eigen::MatrixXd rxhat = so::from_vector({6.0}).get_matrix();
    x2.set_vector({4.0, 5.0, 6.0});
    x2hat = x2.get_matrix();
    REQUIRE(x2hat.rows() == 3);
    REQUIRE(x2hat.cols() == 3);
    CHECK(x2hat(0, 2) == 4.0);
    CHECK(x2hat(1, 2) == 5.0);

    CHECK(x2hat(0, 0) == rxhat(0, 0));
    CHECK(x2hat(0, 1) == rxhat(0, 1));
    CHECK(x2hat(1, 0) == rxhat(1, 0));
    CHECK(x2hat(1, 1) == rxhat(1, 1));
    CHECK(x2hat(2, 0) == 0.0);
    CHECK(x2hat(2, 1) == 0.0);
    CHECK(x2hat(2, 2) == 0.0);
}

TEST_CASE("se operator()", "[domain]")
{
    using namespace Lielab::domain;

    se xzero = se::from_shape(0);
    xzero.set_vector({});

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

    se x1 = se(1);
    x1.set_vector({1.0});

    // In bounds
    CHECK(x1(0) == 1.0);
    CHECK(x1(-1) == 1.0);

    // Out of bounds
    CHECK(std::isnan(x1(-2)));
    CHECK(std::isnan(x1(1)));

    // In bounds
    CHECK(x1(0, 0) == 0.0);
    CHECK(x1(0, 1) == 1.0);
    CHECK(x1(1, 0) == 0.0);
    CHECK(x1(1, 1) == 0.0);
    CHECK(x1(-1, -1) == 0.0);
    CHECK(x1(-1, -2) == 0.0);
    CHECK(x1(-2, -1) == 1.0);
    CHECK(x1(-2, -2) == 0.0);

    // Out of bounds
    CHECK(std::isnan(x1(0, -3)));
    CHECK(std::isnan(x1(-3, 0)));
    CHECK(std::isnan(x1(-3, -3)));
    CHECK(std::isnan(x1(0, 2)));
    CHECK(std::isnan(x1(2, 0)));
    CHECK(std::isnan(x1(2, 2)));

    const Eigen::MatrixXd rxhat = so::from_vector({3.0}).get_matrix();
    se x2 = se(2);
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
    CHECK(x2(0, 0) == rxhat(0, 0));
    CHECK(x2(0, 1) == rxhat(0, 1));
    CHECK(x2(0, 2) == 1.0);
    CHECK(x2(1, 0) == rxhat(1, 0));
    CHECK(x2(1, 1) == rxhat(1, 1));
    CHECK(x2(1, 2) == 2.0);
    CHECK(x2(2, 0) == 0.0);
    CHECK(x2(2, 1) == 0.0);
    CHECK(x2(2, 2) == 0.0);
    CHECK(x2(-1, -1) == 0.0);
    CHECK(x2(-1, -2) == 0.0);
    CHECK(x2(-1, -3) == 0.0);
    CHECK(x2(-2, -1) == 2.0);
    CHECK(x2(-2, -2) == rxhat(1, 1));
    CHECK(x2(-2, -3) == rxhat(1, 0));
    CHECK(x2(-3, -1) == 1.0);
    CHECK(x2(-3, -2) == rxhat(0, 1));
    CHECK(x2(-3, -3) == rxhat(0, 0));

    // Out of bounds
    CHECK(std::isnan(x2(0, -4)));
    CHECK(std::isnan(x2(-4, 0)));
    CHECK(std::isnan(x2(-4, -4)));
    CHECK(std::isnan(x2(0, 3)));
    CHECK(std::isnan(x2(3, 0)));
    CHECK(std::isnan(x2(3, 3)));
}

// TODO: Math ops int

TEST_CASE("se math_ops_double", "[domain]")
{
    using namespace Lielab::domain;

    se x1(2);
    x1.set_vector({1.25, 2.5, 3.75});

    const se x1_lm_2 = 2.0*x1;
    CHECK(x1_lm_2(0) == 2.5);
    CHECK(x1_lm_2(1) == 5.0);
    CHECK(x1_lm_2(2) == 7.5);

    const se x1_rm_2 = x1*2.0;
    CHECK(x1_rm_2(0) == 2.5);
    CHECK(x1_rm_2(1) == 5.0);
    CHECK(x1_rm_2(2) == 7.5);

    x1 *= 2.0;
    CHECK(x1(0) == 2.5);
    CHECK(x1(1) == 5.0);
    CHECK(x1(2) == 7.5);

    x1.set_vector({1.25, 2.5, 3.75});

    const se x1_d_2 = x1/2.0;
    CHECK(x1_d_2(0) == 0.625);
    CHECK(x1_d_2(1) == 1.25);
    CHECK(x1_d_2(2) == 1.875);

    x1 /= 2.0;
    CHECK(x1(0) == 0.625);
    CHECK(x1(1) == 1.25);
    CHECK(x1(2) == 1.875);
}

TEST_CASE("se math_ops_se", "[domain]")
{
    using namespace Lielab::domain;

    se x1(2), x2(2);
    x1.set_vector({1.0, 2.0, 3.0});
    x2.set_vector({1.25, 2.5, 3.75});

    const se x1_add_x2 = x1 + x2;
    CHECK(x1_add_x2(0) == 2.25);
    CHECK(x1_add_x2(1) == 4.5);
    CHECK(x1_add_x2(2) == 6.75);

    x1 += x2;
    CHECK(x1(0) == 2.25);
    CHECK(x1(1) == 4.5);
    CHECK(x1(2) == 6.75);

    x1.set_vector({1.0, 2.0, 3.0});

    const se x1_sub_x2 = x1 - x2;
    CHECK(x1_sub_x2(0) == -0.25);
    CHECK(x1_sub_x2(1) == -0.5);
    CHECK(x1_sub_x2(2) == -0.75);

    x1 -= x2;
    CHECK(x1(0) == -0.25);
    CHECK(x1(1) == -0.5);
    CHECK(x1(2) == -0.75);

    x1.set_vector({1.0, 2.0, 3.0});

    const se x1_unary_sub = (-x1);
    CHECK(x1_unary_sub(0) == -1.0);
    CHECK(x1_unary_sub(1) == -2.0);
    CHECK(x1_unary_sub(2) == -3.0);
}

// TODO: Test projection once it uses so projection
// TEST_CASE("se project", "[domain]")
// {
//     using namespace Lielab::domain;

//     const Eigen::MatrixXd rand_2_2 = Eigen::MatrixXd::Random(2, 2);
//     const Eigen::MatrixXd proj_2_2 = rn::project(rand_2_2);

//     REQUIRE(proj_2_2.rows() == 2);
//     REQUIRE(proj_2_2.cols() == 2);
//     CHECK(proj_2_2(0, 0) == 0.0);
//     CHECK(proj_2_2(0, 1) == rand_2_2(0, 1));
//     CHECK(proj_2_2(1, 0) == 0.0);
//     CHECK(proj_2_2(1, 1) == 0.0);

//     const Eigen::MatrixXd rand_3_3 = Eigen::MatrixXd::Random(3, 3);
//     const Eigen::MatrixXd proj_3_3 = rn::project(rand_3_3);

//     REQUIRE(proj_3_3.rows() == 3);
//     REQUIRE(proj_3_3.cols() == 3);
//     CHECK(proj_3_3(0, 0) == 0.0);
//     CHECK(proj_3_3(0, 1) == 0.0);
//     CHECK(proj_3_3(0, 2) == rand_3_3(0, 2));
//     CHECK(proj_3_3(1, 0) == 0.0);
//     CHECK(proj_3_3(1, 1) == 0.0);
//     CHECK(proj_3_3(1, 2) == rand_3_3(1, 2));
//     CHECK(proj_3_3(2, 0) == 0.0);
//     CHECK(proj_3_3(2, 1) == 0.0);
//     CHECK(proj_3_3(2, 2) == 0.0);

//     const Eigen::MatrixXd rand_2_3 = Eigen::MatrixXd::Random(2, 3);
//     const Eigen::MatrixXd proj_2_3 = rn::project(rand_2_3);

//     REQUIRE(proj_2_3.rows() == 2);
//     REQUIRE(proj_2_3.cols() == 2);
//     CHECK(proj_2_3(0, 0) == 0.0);
//     CHECK(proj_2_3(0, 1) == rand_2_3(0, 1));
//     CHECK(proj_2_3(1, 0) == 0.0);
//     CHECK(proj_2_3(1, 1) == 0.0);

//     const Eigen::MatrixXd rand_3_2 = Eigen::MatrixXd::Random(3, 2);
//     const Eigen::MatrixXd proj_3_2 = rn::project(rand_3_2);

//     REQUIRE(proj_3_2.rows() == 2);
//     REQUIRE(proj_3_2.cols() == 2);
//     CHECK(proj_3_2(0, 0) == 0.0);
//     CHECK(proj_3_2(0, 1) == rand_3_2(0, 1));
//     CHECK(proj_3_2(1, 0) == 0.0);
//     CHECK(proj_3_2(1, 1) == 0.0);
// }

TEST_CASE("se get/from_vector", "[domain]")
{
    /*!
    * Tests the get/from_vector operation.
    */

    using namespace Lielab::domain;

    const se x0 = se::from_vector({});
    const Eigen::VectorXd x0bar = x0.get_vector();

    CHECK(x0.get_shape() == 1);
    CHECK(x0bar.size() == 0);

    const se x1 = se::from_vector({1.0});
    const Eigen::VectorXd x1bar = x1.get_vector();

    CHECK(x1.get_shape() == 2);
    REQUIRE(x1bar.size() == 1);
    CHECK(x1bar(0) == 1.0);

    const se x2 = se::from_vector({1.0, 2.0});
    const Eigen::VectorXd x2bar = x2.get_vector();

    CHECK(x2.get_shape() == 3);
    REQUIRE(x2bar.size() == 3);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 0.0);

    const se x3 = se::from_vector({1.0, 2.0, 3.0});
    const Eigen::VectorXd x3bar = x3.get_vector();

    CHECK(x3.get_shape() == 3);
    REQUIRE(x3bar.size() == 3);
    CHECK(x3bar(0) == 1.0);
    CHECK(x3bar(1) == 2.0);
    CHECK(x3bar(2) == 3.0);
}

// TEST_CASE("se4", "[domain]")
// {
//     /*!
//     * Tests the se algebra with se(4).
//     */

//     Lielab::domain::se x = Lielab::domain::se::basis(0, 4);
//     Lielab::domain::se y = Lielab::domain::se::basis(1, 4);
//     Lielab::domain::se z = Lielab::domain::se::basis(2, 4);
//     Lielab::domain::se u = Lielab::domain::se::basis(3, 4);
//     Lielab::domain::se v = Lielab::domain::se::basis(4, 4);
//     Lielab::domain::se w = Lielab::domain::se::basis(5, 4);
//     Lielab::domain::se zero = x*0;

//     assert_domain(Lielab::functions::commutator(x, y), zero);
//     assert_domain(Lielab::functions::commutator(y, z), zero);
//     assert_domain(Lielab::functions::commutator(z, x), zero);
//     assert_domain(Lielab::functions::commutator(y, x), zero);
//     assert_domain(Lielab::functions::commutator(z, y), zero);
//     assert_domain(Lielab::functions::commutator(x, z), zero);
//     assert_domain(Lielab::functions::commutator(x, v), z);
//     assert_domain(Lielab::functions::commutator(y, w), x);
//     assert_domain(Lielab::functions::commutator(z, u), y);
//     assert_domain(Lielab::functions::commutator(y, u), -z);
//     assert_domain(Lielab::functions::commutator(z, v), -x);
//     assert_domain(Lielab::functions::commutator(x, w), -y);
//     assert_domain(Lielab::functions::commutator(u, v), w);
//     assert_domain(Lielab::functions::commutator(v, w), u);
//     assert_domain(Lielab::functions::commutator(w, u), v);
//     assert_domain(Lielab::functions::commutator(v, u), -w);
//     assert_domain(Lielab::functions::commutator(w, v), -u);
//     assert_domain(Lielab::functions::commutator(u, w), -v);
// }
