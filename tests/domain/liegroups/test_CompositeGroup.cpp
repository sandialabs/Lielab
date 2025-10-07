#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

#include <iostream>

Lielab::domain::CompositeGroup _make_cgroup()
{
    using namespace Lielab::domain;

    Lielab::domain::CN yCN1 = Lielab::domain::CN(2);
    yCN1.unserialize({1.0, 2.0, 3.0, 4.0});
    Lielab::domain::GLC yGLC1 = Lielab::domain::GLC(2);
    yGLC1.unserialize({5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0});
    Lielab::domain::GLR yGLR1 = Lielab::domain::GLR(2);
    yGLR1.unserialize({13.0, 14.0, 15.0, 16.0});
    Lielab::domain::RN yRN1 = Lielab::domain::RN(2);
    yRN1.unserialize({17.0, 18.0});
    Lielab::domain::SE ySE1 = Lielab::domain::SE(2);
    ySE1.unserialize({19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0});
    Lielab::domain::SO ySO1 = Lielab::domain::SO(2);
    ySO1.unserialize({28.0, 29.0, 30.0, 31.0});
    Lielab::domain::SP ySP1 = Lielab::domain::SP(2);
    ySP1.unserialize({32.0, 33.0, 34.0, 35.0});
    Lielab::domain::SU ySU1 = Lielab::domain::SU(2);
    ySU1.unserialize({36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0});

    return CompositeGroup({yCN1, yGLC1, yGLR1, yRN1, ySE1, ySO1, ySP1, ySU1});
}

TEST_CASE("CompositeGroup to_string", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeGroup y1 = _make_cgroup();

    CHECK(y1.to_string() == "C^2 x GL(2, C) x GL(2, R) x R^2 x SE(2) x SO(2) x SP(2, R) x SU(2)");
}

TEST_CASE("CompositeGroup main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeGroup xblank = CompositeGroup();
    CHECK(xblank.get_dimension() == 0);

    const CompositeGroup x0 = CompositeGroup(0);
    CHECK(x0.get_dimension() == 0);
    const CompositeGroup x1 = CompositeGroup(1);
    CHECK(x1.get_dimension() == 2);
    const CompositeGroup x10 = CompositeGroup(10);
    CHECK(x10.get_dimension() == 200);
}

TEST_CASE("CompositeGroup from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeGroup x0 = CompositeGroup::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXcd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const CompositeGroup x1 = CompositeGroup::from_shape(1);
    CHECK(x1.get_dimension() == 2);
    const Eigen::MatrixXcd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const CompositeGroup x2 = CompositeGroup::from_shape(2);
    CHECK(x2.get_dimension() == 8);
    const Eigen::MatrixXcd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("CompositeGroup get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    CompositeGroup zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

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

    CompositeGroup y1 = _make_cgroup();
    
    const std::vector<size_t> dims = y1.get_dimensions();
    REQUIRE(dims.size() == 8);
    CHECK(dims[0] == 4);
    CHECK(dims[1] == 8);
    CHECK(dims[2] == 4);
    CHECK(dims[3] == 2);
    CHECK(dims[4] == 3);
    CHECK(dims[5] == 1);
    CHECK(dims[6] == 3);
    CHECK(dims[7] == 3);

    CHECK(y1.get_dimension() == 28);
}

TEST_CASE("CompositeGroup get_size", "[domain]")
{
    using namespace Lielab::domain;

    CompositeGroup zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(zero.get_size() == 0);
    CHECK(one.get_size() == 2);
    CHECK(two.get_size() == 8);
    CHECK(three.get_size() == 18);
    CHECK(four.get_size() == 32);
    CHECK(five.get_size() == 50);
    CHECK(six.get_size() == 72);
    CHECK(seven.get_size() == 98);
    CHECK(eight.get_dimension() == 128);

    CompositeGroup y1 = _make_cgroup();
    
    const std::vector<size_t> sizes = y1.get_sizes();
    REQUIRE(sizes.size() == 8);
    CHECK(sizes[0] == 4);
    CHECK(sizes[1] == 8);
    CHECK(sizes[2] == 4);
    CHECK(sizes[3] == 2);
    CHECK(sizes[4] == 9);
    CHECK(sizes[5] == 4);
    CHECK(sizes[6] == 4);
    CHECK(sizes[7] == 8);

    CHECK(y1.get_size() == 43);
}

TEST_CASE("CompositeGroup serialize/unserialize", "[domain]")
{
    /*!
    * Tests the serialize/unserialize operation.
    */

    using namespace Lielab::domain;

    CompositeGroup y1 = _make_cgroup();
    y1.unserialize({43.0, 42.0, 41.0, 40.0, 39.0, 38.0, 37.0, 36.0, 35.0, 34.0, 33.0, 32.0, 31.0, 30.0,
                    29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0, 20.0, 19.0, 18.0, 17.0, 16.0,
                    15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0});
    
    const Eigen::VectorXd y1bar = y1.serialize();
    REQUIRE(y1bar.size() == 43);
    CHECK(y1bar(0) == 43.0);
    CHECK(y1bar(1) == 42.0);
    CHECK(y1bar(2) == 41.0);
    CHECK(y1bar(3) == 40.0);
    CHECK(y1bar(4) == 39.0);
    CHECK(y1bar(5) == 38.0);
    CHECK(y1bar(6) == 37.0);
    CHECK(y1bar(7) == 36.0);
    CHECK(y1bar(8) == 35.0);
    CHECK(y1bar(9) == 34.0);
    CHECK(y1bar(10) == 33.0);
    CHECK(y1bar(11) == 32.0);
    CHECK(y1bar(12) == 31.0);
    CHECK(y1bar(13) == 30.0);
    CHECK(y1bar(14) == 29.0);
    CHECK(y1bar(15) == 28.0);
    CHECK(y1bar(16) == 27.0);
    CHECK(y1bar(17) == 26.0);
    CHECK(y1bar(18) == 25.0);
    CHECK(y1bar(19) == 24.0);
    CHECK(y1bar(20) == 23.0);
    CHECK(y1bar(21) == 22.0);
    CHECK(y1bar(22) == 21.0);
    CHECK(y1bar(23) == 20.0);
    CHECK(y1bar(24) == 19.0);
    CHECK(y1bar(25) == 18.0);
    CHECK(y1bar(26) == 17.0);
    CHECK(y1bar(27) == 16.0);
    CHECK(y1bar(28) == 15.0);
    CHECK(y1bar(29) == 14.0);
    CHECK(y1bar(30) == 13.0);
    CHECK(y1bar(31) == 12.0);
    CHECK(y1bar(32) == 11.0);
    CHECK(y1bar(33) == 10.0);
    CHECK(y1bar(34) == 9.0);
    CHECK(y1bar(35) == 8.0);
    CHECK(y1bar(36) == 7.0);
    CHECK(y1bar(37) == 6.0);
    CHECK(y1bar(38) == 5.0);
    CHECK(y1bar(39) == 4.0);
    CHECK(y1bar(40) == 3.0);
    CHECK(y1bar(41) == 2.0);
    CHECK(y1bar(42) == 1.0);
}

TEST_CASE("CompositeGroup get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    const CompositeGroup y1 = _make_cgroup();

    const Eigen::MatrixXcd y1hat = y1.get_matrix();

    REQUIRE(y1hat.rows() == 19);
    REQUIRE(y1hat.cols() == 19);
    
    // TODO: Check the 0's
    // TODO: Check each submatrix

}

TEST_CASE("CompositeGroup operator()", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeGroup y1 = _make_cgroup();

    // In bounds CN component
    CHECK(y1(0, 0) == std::complex<double>(1.0, 0.0));
    CHECK(y1(0, 1) == std::complex<double>(0.0, 0.0));
    CHECK(y1(0, 2) == std::complex<double>(1.0, 2.0));
    CHECK(y1(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(y1(1, 1) == std::complex<double>(1.0, 0.0));
    CHECK(y1(1, 2) == std::complex<double>(3.0, 4.0));
    CHECK(y1(2, 0) == std::complex<double>(0.0, 0.0));
    CHECK(y1(2, 1) == std::complex<double>(0.0, 0.0));
    CHECK(y1(2, 2) == std::complex<double>(1.0, 0.0));
    CHECK(y1(-19, -19) == std::complex<double>(1.0, 0.0));
    CHECK(y1(-19, -18) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-19, -17) == std::complex<double>(1.0, 2.0));
    CHECK(y1(-18, -19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-18, -18) == std::complex<double>(1.0, 0.0));
    CHECK(y1(-18, -17) == std::complex<double>(3.0, 4.0));
    CHECK(y1(-17, -19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-17, -18) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-17, -17) == std::complex<double>(1.0, 0.0));

    // In bounds GLC component
    CHECK(y1(3, 3) == std::complex<double>(5.0, 6.0));
    CHECK(y1(3, 4) == std::complex<double>(7.0, 8.0));
    CHECK(y1(4, 3) == std::complex<double>(9.0, 10.0));
    CHECK(y1(4, 4) == std::complex<double>(11.0, 12.0));
    CHECK(y1(-16, -16) == std::complex<double>(5.0, 6.0));
    CHECK(y1(-16, -15) == std::complex<double>(7.0, 8.0));
    CHECK(y1(-15, -16) == std::complex<double>(9.0, 10.0));
    CHECK(y1(-15, -15) == std::complex<double>(11.0, 12.0));

    // In bounds GLR component
    CHECK(y1(5, 5) == std::complex<double>(13.0, 0.0));
    CHECK(y1(5, 6) == std::complex<double>(14.0, 0.0));
    CHECK(y1(6, 5) == std::complex<double>(15.0, 0.0));
    CHECK(y1(6, 6) == std::complex<double>(16.0, 0.0));
    CHECK(y1(-14, -14) == std::complex<double>(13.0, 0.0));
    CHECK(y1(-14, -13) == std::complex<double>(14.0, 0.0));
    CHECK(y1(-13, -14) == std::complex<double>(15.0, 0.0));
    CHECK(y1(-13, -13) == std::complex<double>(16.0, 0.0));

    // In bounds RN component
    CHECK(y1(7, 7) == std::complex<double>(1.0, 0.0));
    CHECK(y1(7, 8) == std::complex<double>(0.0, 0.0));
    CHECK(y1(7, 9) == std::complex<double>(17.0, 0.0));
    CHECK(y1(8, 7) == std::complex<double>(0.0, 0.0));
    CHECK(y1(8, 8) == std::complex<double>(1.0, 0.0));
    CHECK(y1(8, 9) == std::complex<double>(18.0, 0.0));
    CHECK(y1(9, 7) == std::complex<double>(0.0, 0.0));
    CHECK(y1(9, 8) == std::complex<double>(0.0, 0.0));
    CHECK(y1(9, 9) == std::complex<double>(1.0, 0.0));
    CHECK(y1(-12, -12) == std::complex<double>(1.0, 0.0));
    CHECK(y1(-12, -11) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-12, -10) == std::complex<double>(17.0, 0.0));
    CHECK(y1(-11, -12) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-11, -11) == std::complex<double>(1.0, 0.0));
    CHECK(y1(-11, -10) == std::complex<double>(18.0, 0.0));
    CHECK(y1(-10, -12) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-10, -11) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-10, -10) == std::complex<double>(1.0, 0.0));

    // In bounds SE component
    CHECK(y1(10, 10) == std::complex<double>(19.0, 0.0));
    CHECK(y1(10, 11) == std::complex<double>(20.0, 0.0));
    CHECK(y1(10, 12) == std::complex<double>(21.0, 0.0));
    CHECK(y1(11, 10) == std::complex<double>(22.0, 0.0));
    CHECK(y1(11, 11) == std::complex<double>(23.0, 0.0));
    CHECK(y1(11, 12) == std::complex<double>(24.0, 0.0));
    CHECK(y1(12, 10) == std::complex<double>(25.0, 0.0));
    CHECK(y1(12, 11) == std::complex<double>(26.0, 0.0));
    CHECK(y1(12, 12) == std::complex<double>(27.0, 0.0));
    CHECK(y1(-9, -9) == std::complex<double>(19.0, 0.0));
    CHECK(y1(-9, -8) == std::complex<double>(20.0, 0.0));
    CHECK(y1(-9, -7) == std::complex<double>(21.0, 0.0));
    CHECK(y1(-8, -9) == std::complex<double>(22.0, 0.0));
    CHECK(y1(-8, -8) == std::complex<double>(23.0, 0.0));
    CHECK(y1(-8, -7) == std::complex<double>(24.0, 0.0));
    CHECK(y1(-7, -9) == std::complex<double>(25.0, 0.0));
    CHECK(y1(-7, -8) == std::complex<double>(26.0, 0.0));
    CHECK(y1(-7, -7) == std::complex<double>(27.0, 0.0));

    // In bounds SO component
    CHECK(y1(13, 13) == std::complex<double>(28.0, 0.0));
    CHECK(y1(13, 14) == std::complex<double>(29.0, 0.0));
    CHECK(y1(14, 13) == std::complex<double>(30.0, 0.0));
    CHECK(y1(14, 14) == std::complex<double>(31.0, 0.0));
    CHECK(y1(-6, -6) == std::complex<double>(28.0, 0.0));
    CHECK(y1(-6, -5) == std::complex<double>(29.0, 0.0));
    CHECK(y1(-5, -6) == std::complex<double>(30.0, 0.0));
    CHECK(y1(-5, -5) == std::complex<double>(31.0, 0.0));

    // In bounds SP component
    CHECK(y1(15, 15) == std::complex<double>(32.0, 0.0));
    CHECK(y1(15, 16) == std::complex<double>(33.0, 0.0));
    CHECK(y1(16, 15) == std::complex<double>(34.0, 0.0));
    CHECK(y1(16, 16) == std::complex<double>(35.0, 0.0));
    CHECK(y1(-4, -4) == std::complex<double>(32.0, 0.0));
    CHECK(y1(-4, -3) == std::complex<double>(33.0, 0.0));
    CHECK(y1(-3, -4) == std::complex<double>(34.0, 0.0));
    CHECK(y1(-3, -3) == std::complex<double>(35.0, 0.0));

    // In bounds SU component
    CHECK(y1(17, 17) == std::complex<double>(36.0, 37.0));
    CHECK(y1(17, 18) == std::complex<double>(38.0, 39.0));
    CHECK(y1(18, 17) == std::complex<double>(40.0, 41.0));
    CHECK(y1(18, 18) == std::complex<double>(42.0, 43.0));
    CHECK(y1(-2, -2) == std::complex<double>(36.0, 37.0));
    CHECK(y1(-2, -1) == std::complex<double>(38.0, 39.0));
    CHECK(y1(-1, -2) == std::complex<double>(40.0, 41.0));
    CHECK(y1(-1, -1) == std::complex<double>(42.0, 43.0));

    // In bounds sparse components
    CHECK(y1(0, 3) == std::complex<double>(0.0, 0.0));
    CHECK(y1(0, 4) == std::complex<double>(0.0, 0.0));
    CHECK(y1(1, 3) == std::complex<double>(0.0, 0.0));
    CHECK(y1(1, 4) == std::complex<double>(0.0, 0.0));
    CHECK(y1(2, 3) == std::complex<double>(0.0, 0.0));
    CHECK(y1(2, 4) == std::complex<double>(0.0, 0.0));
    CHECK(y1(3, 0) == std::complex<double>(0.0, 0.0));
    CHECK(y1(4, 0) == std::complex<double>(0.0, 0.0));
    CHECK(y1(3, 1) == std::complex<double>(0.0, 0.0));
    CHECK(y1(4, 1) == std::complex<double>(0.0, 0.0));
    CHECK(y1(3, 2) == std::complex<double>(0.0, 0.0));
    CHECK(y1(4, 2) == std::complex<double>(0.0, 0.0));
    // Probably don't need to check the rest

    // Out of bounds
    CHECK(std::isnan(y1(0, -20).real()));
    CHECK(std::isnan(y1(0, -20).imag()));
    CHECK(std::isnan(y1(-20, 0).real()));
    CHECK(std::isnan(y1(-20, 0).imag()));
    CHECK(std::isnan(y1(-20, -20).real()));
    CHECK(std::isnan(y1(-20, -20).imag()));
    CHECK(std::isnan(y1(0, 19).real()));
    CHECK(std::isnan(y1(0, 19).imag()));
    CHECK(std::isnan(y1(19, 0).real()));
    CHECK(std::isnan(y1(19, 0).imag()));
    CHECK(std::isnan(y1(19, 19).real()));
    CHECK(std::isnan(y1(19, 19).imag()));
}

TEST_CASE("CompositeGroup math_ops_CompositeGroup", "[domain]")
{
    using namespace Lielab::domain;

    CompositeGroup y1 = _make_cgroup();
    CompositeGroup y2 = _make_cgroup();

    const CompositeGroup y1_prod_y2 = y1*y2;
    REQUIRE(y1_prod_y2.space.size() == y1.space.size());

    const Eigen::VectorXd y1_prod_y2_0bar = std::get<CN>(y1_prod_y2.space[0]).serialize();
    REQUIRE(y1_prod_y2_0bar.size() == 4);
    CHECK(y1_prod_y2_0bar(0) == 2.0);
    CHECK(y1_prod_y2_0bar(1) == 4.0);
    CHECK(y1_prod_y2_0bar(2) == 6.0);
    CHECK(y1_prod_y2_0bar(3) == 8.0);

    const Eigen::VectorXd y1_prod_y2_1bar = std::get<GLC>(y1_prod_y2.space[1]).serialize();
    REQUIRE(y1_prod_y2_1bar.size() == 8);
    CHECK(y1_prod_y2_1bar(0) == -28.0);
    CHECK(y1_prod_y2_1bar(1) == 202.0);
    CHECK(y1_prod_y2_1bar(2) == -32.0);
    CHECK(y1_prod_y2_1bar(3) == 254.0);
    CHECK(y1_prod_y2_1bar(4) == -36.0);
    CHECK(y1_prod_y2_1bar(5) == 322.0);
    CHECK(y1_prod_y2_1bar(6) == -40.0);
    CHECK(y1_prod_y2_1bar(7) == 406.0);

    const Eigen::VectorXd y1_prod_y2_2bar = std::get<GLR>(y1_prod_y2.space[2]).serialize();
    REQUIRE(y1_prod_y2_2bar.size() == 4);
    CHECK(y1_prod_y2_2bar(0) == 379.0);
    CHECK(y1_prod_y2_2bar(1) == 406.0);
    CHECK(y1_prod_y2_2bar(2) == 435.0);
    CHECK(y1_prod_y2_2bar(3) == 466.0);

    const Eigen::VectorXd y1_prod_y2_3bar = std::get<RN>(y1_prod_y2.space[3]).serialize();
    REQUIRE(y1_prod_y2_3bar.size() == 2);
    CHECK(y1_prod_y2_3bar(0) == 34.0);
    CHECK(y1_prod_y2_3bar(1) == 36.0);

    const Eigen::VectorXd y1_prod_y2_4bar = std::get<SE>(y1_prod_y2.space[4]).serialize();
    REQUIRE(y1_prod_y2_4bar.size() == 9);
    CHECK(y1_prod_y2_4bar(0) == 1326.0);
    CHECK(y1_prod_y2_4bar(1) == 1386.0);
    CHECK(y1_prod_y2_4bar(2) == 1446.0);
    CHECK(y1_prod_y2_4bar(3) == 1524.0);
    CHECK(y1_prod_y2_4bar(4) == 1593.0);
    CHECK(y1_prod_y2_4bar(5) == 1662.0);
    CHECK(y1_prod_y2_4bar(6) == 1722.0);
    CHECK(y1_prod_y2_4bar(7) == 1800.0);
    CHECK(y1_prod_y2_4bar(8) == 1878.0);

    const Eigen::VectorXd y1_prod_y2_5bar = std::get<SO>(y1_prod_y2.space[5]).serialize();
    REQUIRE(y1_prod_y2_5bar.size() == 4);
    CHECK(y1_prod_y2_5bar(0) == 1654.0);
    CHECK(y1_prod_y2_5bar(1) == 1711.0);
    CHECK(y1_prod_y2_5bar(2) == 1770.0);
    CHECK(y1_prod_y2_5bar(3) == 1831.0);

    const Eigen::VectorXd y1_prod_y2_6bar = std::get<SP>(y1_prod_y2.space[6]).serialize();
    REQUIRE(y1_prod_y2_6bar.size() == 4);
    CHECK(y1_prod_y2_6bar(0) == 2146.0);
    CHECK(y1_prod_y2_6bar(1) == 2211.0);
    CHECK(y1_prod_y2_6bar(2) == 2278.0);
    CHECK(y1_prod_y2_6bar(3) == 2347.0);

    const Eigen::VectorXd y1_prod_y2_7bar = std::get<SU>(y1_prod_y2.space[7]).serialize();
    REQUIRE(y1_prod_y2_7bar.size() == 8);
    CHECK(y1_prod_y2_7bar(0) == -152.0);
    CHECK(y1_prod_y2_7bar(1) == 5782.0);
    CHECK(y1_prod_y2_7bar(2) == -156.0);
    CHECK(y1_prod_y2_7bar(3) == 6082.0);
    CHECK(y1_prod_y2_7bar(4) == -160.0);
    CHECK(y1_prod_y2_7bar(5) == 6398.0);
    CHECK(y1_prod_y2_7bar(6) == -164.0);
    CHECK(y1_prod_y2_7bar(7) == 6730.0);

    y1 *= y2;
    REQUIRE(y1.space.size() == y2.space.size());

    const Eigen::VectorXd y1_iprod_y2_0bar = std::get<CN>(y1.space[0]).serialize();
    REQUIRE(y1_iprod_y2_0bar.size() == 4);
    CHECK(y1_iprod_y2_0bar(0) == 2.0);
    CHECK(y1_iprod_y2_0bar(1) == 4.0);
    CHECK(y1_iprod_y2_0bar(2) == 6.0);
    CHECK(y1_iprod_y2_0bar(3) == 8.0);

    const Eigen::VectorXd y1_iprod_y2_1bar = std::get<GLC>(y1.space[1]).serialize();
    REQUIRE(y1_iprod_y2_1bar.size() == 8);
    CHECK(y1_iprod_y2_1bar(0) == -28.0);
    CHECK(y1_iprod_y2_1bar(1) == 202.0);
    CHECK(y1_iprod_y2_1bar(2) == -32.0);
    CHECK(y1_iprod_y2_1bar(3) == 254.0);
    CHECK(y1_iprod_y2_1bar(4) == -36.0);
    CHECK(y1_iprod_y2_1bar(5) == 322.0);
    CHECK(y1_iprod_y2_1bar(6) == -40.0);
    CHECK(y1_iprod_y2_1bar(7) == 406.0);

    const Eigen::VectorXd y1_iprod_y2_2bar = std::get<GLR>(y1.space[2]).serialize();
    REQUIRE(y1_iprod_y2_2bar.size() == 4);
    CHECK(y1_iprod_y2_2bar(0) == 379.0);
    CHECK(y1_iprod_y2_2bar(1) == 406.0);
    CHECK(y1_iprod_y2_2bar(2) == 435.0);
    CHECK(y1_iprod_y2_2bar(3) == 466.0);

    const Eigen::VectorXd y1_iprod_y2_3bar = std::get<RN>(y1.space[3]).serialize();
    REQUIRE(y1_iprod_y2_3bar.size() == 2);
    CHECK(y1_iprod_y2_3bar(0) == 34.0);
    CHECK(y1_iprod_y2_3bar(1) == 36.0);

    const Eigen::VectorXd y1_iprod_y2_4bar = std::get<SE>(y1.space[4]).serialize();
    REQUIRE(y1_iprod_y2_4bar.size() == 9);
    CHECK(y1_iprod_y2_4bar(0) == 1326.0);
    CHECK(y1_iprod_y2_4bar(1) == 1386.0);
    CHECK(y1_iprod_y2_4bar(2) == 1446.0);
    CHECK(y1_iprod_y2_4bar(3) == 1524.0);
    CHECK(y1_iprod_y2_4bar(4) == 1593.0);
    CHECK(y1_iprod_y2_4bar(5) == 1662.0);
    CHECK(y1_iprod_y2_4bar(6) == 1722.0);
    CHECK(y1_iprod_y2_4bar(7) == 1800.0);
    CHECK(y1_iprod_y2_4bar(8) == 1878.0);

    const Eigen::VectorXd y1_iprod_y2_5bar = std::get<SO>(y1.space[5]).serialize();
    REQUIRE(y1_iprod_y2_5bar.size() == 4);
    CHECK(y1_iprod_y2_5bar(0) == 1654.0);
    CHECK(y1_iprod_y2_5bar(1) == 1711.0);
    CHECK(y1_iprod_y2_5bar(2) == 1770.0);
    CHECK(y1_iprod_y2_5bar(3) == 1831.0);

    const Eigen::VectorXd y1_iprod_y2_6bar = std::get<SP>(y1.space[6]).serialize();
    REQUIRE(y1_iprod_y2_6bar.size() == 4);
    CHECK(y1_iprod_y2_6bar(0) == 2146.0);
    CHECK(y1_iprod_y2_6bar(1) == 2211.0);
    CHECK(y1_iprod_y2_6bar(2) == 2278.0);
    CHECK(y1_iprod_y2_6bar(3) == 2347.0);

    const Eigen::VectorXd y1_iprod_y2_7bar = std::get<SU>(y1.space[7]).serialize();
    REQUIRE(y1_iprod_y2_7bar.size() == 8);
    CHECK(y1_iprod_y2_7bar(0) == -152.0);
    CHECK(y1_iprod_y2_7bar(1) == 5782.0);
    CHECK(y1_iprod_y2_7bar(2) == -156.0);
    CHECK(y1_iprod_y2_7bar(3) == 6082.0);
    CHECK(y1_iprod_y2_7bar(4) == -160.0);
    CHECK(y1_iprod_y2_7bar(5) == 6398.0);
    CHECK(y1_iprod_y2_7bar(6) == -164.0);
    CHECK(y1_iprod_y2_7bar(7) == 6730.0);

    y1 = _make_cgroup();

    const CompositeGroup y1_inv = y1.inverse();
    REQUIRE(y1_inv.space.size() == y1.space.size());

    const Eigen::VectorXd y1_inv_0bar = std::get<CN>(y1_inv.space[0]).serialize();
    REQUIRE(y1_inv_0bar.size() == 4);
    CHECK(y1_inv_0bar(0) == -1.0);
    CHECK(y1_inv_0bar(1) == -2.0);
    CHECK(y1_inv_0bar(2) == -3.0);
    CHECK(y1_inv_0bar(3) == -4.0);

    const Eigen::VectorXd y1_inv_1bar = std::get<GLC>(y1_inv.space[1]).serialize();
    REQUIRE(y1_inv_1bar.size() == 8);
    CHECK_THAT(y1_inv_1bar(0), Catch::Matchers::WithinAbs(-0.75, 1e-14));
    CHECK_THAT(y1_inv_1bar(1), Catch::Matchers::WithinAbs(0.6875, 1e-14));
    CHECK_THAT(y1_inv_1bar(2), Catch::Matchers::WithinAbs(0.5, 1e-14));
    CHECK_THAT(y1_inv_1bar(3), Catch::Matchers::WithinAbs(-0.4375, 1e-14));
    CHECK_THAT(y1_inv_1bar(4), Catch::Matchers::WithinAbs(0.625, 1e-14));
    CHECK_THAT(y1_inv_1bar(5), Catch::Matchers::WithinAbs(-0.5625, 1e-14));
    CHECK_THAT(y1_inv_1bar(6), Catch::Matchers::WithinAbs(-0.375, 1e-14));
    CHECK_THAT(y1_inv_1bar(7), Catch::Matchers::WithinAbs(0.3125, 1e-14));

    const Eigen::VectorXd y1_inv_2bar = std::get<GLR>(y1_inv.space[2]).serialize();
    REQUIRE(y1_inv_2bar.size() == 4);
    CHECK_THAT(y1_inv_2bar(0), Catch::Matchers::WithinAbs(-8.0, 1e-13));
    CHECK_THAT(y1_inv_2bar(1), Catch::Matchers::WithinAbs(7.0, 1e-13));
    CHECK_THAT(y1_inv_2bar(2), Catch::Matchers::WithinAbs(7.5, 1e-13));
    CHECK_THAT(y1_inv_2bar(3), Catch::Matchers::WithinAbs(-6.5, 1e-13));

    const Eigen::VectorXd y1_inv_3bar = std::get<RN>(y1_inv.space[3]).serialize();
    REQUIRE(y1_inv_3bar.size() == 2);
    CHECK(y1_inv_3bar(0) == -17.0);
    CHECK(y1_inv_3bar(1) == -18.0);

    const Eigen::VectorXd y1_inv_4bar = std::get<SE>(y1_inv.space[4]).serialize();
    REQUIRE(y1_inv_4bar.size() == 9);
    // TODO: Change default values to something more reasonable
    // CHECK(y1_inv_4bar(0) == 1326.0);
    // CHECK(y1_inv_4bar(1) == 1386.0);
    // CHECK(y1_inv_4bar(2) == 1446.0);
    // CHECK(y1_inv_4bar(3) == 1524.0);
    // CHECK(y1_inv_4bar(4) == 1593.0);
    // CHECK(y1_inv_4bar(5) == 1662.0);
    // CHECK(y1_inv_4bar(6) == 1722.0);
    // CHECK(y1_inv_4bar(7) == 1800.0);
    // CHECK(y1_inv_4bar(8) == 1878.0);

    const Eigen::VectorXd y1_inv_5bar = std::get<SO>(y1_inv.space[5]).serialize();
    REQUIRE(y1_inv_5bar.size() == 4);
    CHECK(y1_inv_5bar(0) == 28.0);
    CHECK(y1_inv_5bar(1) == 30.0);
    CHECK(y1_inv_5bar(2) == 29.0);
    CHECK(y1_inv_5bar(3) == 31.0);
    
    const Eigen::VectorXd y1_inv_6bar = std::get<SP>(y1_inv.space[6]).serialize();
    REQUIRE(y1_inv_6bar.size() == 4);
    CHECK_THAT(y1_inv_6bar(0), Catch::Matchers::WithinAbs(-17.5, 1e-12));
    CHECK_THAT(y1_inv_6bar(1), Catch::Matchers::WithinAbs(16.5, 1e-12));
    CHECK_THAT(y1_inv_6bar(2), Catch::Matchers::WithinAbs(17.0, 1e-12));
    CHECK_THAT(y1_inv_6bar(3), Catch::Matchers::WithinAbs(-16.0, 1e-12));

    const Eigen::VectorXd y1_inv_7bar = std::get<SU>(y1_inv.space[7]).serialize();
    REQUIRE(y1_inv_7bar.size() == 8);
    CHECK_THAT(y1_inv_7bar(0), Catch::Matchers::WithinAbs(-2.6875, 1e-13));
    CHECK_THAT(y1_inv_7bar(1), Catch::Matchers::WithinAbs(2.625, 1e-13));
    CHECK_THAT(y1_inv_7bar(2), Catch::Matchers::WithinAbs(2.4375, 1e-13));
    CHECK_THAT(y1_inv_7bar(3), Catch::Matchers::WithinAbs(-2.375, 1e-13));
    CHECK_THAT(y1_inv_7bar(4), Catch::Matchers::WithinAbs(2.5625, 1e-13));
    CHECK_THAT(y1_inv_7bar(5), Catch::Matchers::WithinAbs(-2.5, 1e-13));
    CHECK_THAT(y1_inv_7bar(6), Catch::Matchers::WithinAbs(-2.3125, 1e-13));
    CHECK_THAT(y1_inv_7bar(7), Catch::Matchers::WithinAbs(2.25, 1e-13));
}

TEST_CASE("CompositeGroup operator[]", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeGroup x1 = _make_cgroup();
    
    // It would be nice if these could be called without std::get<> like:
    // const CN x10 = x1[0];
    const CN x10 = std::get<CN>(x1[0]);
    const GLC x11 = std::get<GLC>(x1[1]);
    const GLR x12 = std::get<GLR>(x1[2]);
    const RN x13 = std::get<RN>(x1[3]);
    const SE x14 = std::get<SE>(x1[4]);
    const SO x15 = std::get<SO>(x1[5]);
    const SP x16 = std::get<SP>(x1[6]);
    const SU x17 = std::get<SU>(x1[7]);
    const CN x1m8 = std::get<CN>(x1[-8]);
    const GLC x1m7 = std::get<GLC>(x1[-7]);
    const GLR x1m6 = std::get<GLR>(x1[-6]);
    const RN x1m5 = std::get<RN>(x1[-5]);
    const SE x1m4 = std::get<SE>(x1[-4]);
    const SO x1m3 = std::get<SO>(x1[-3]);
    const SP x1m2 = std::get<SP>(x1[-2]);
    const SU x1m1 = std::get<SU>(x1[-1]);

    CHECK(x10.to_string() == "C^2");
    CHECK(x11.to_string() == "GL(2, C)");
    CHECK(x12.to_string() == "GL(2, R)");
    CHECK(x13.to_string() == "R^2");
    CHECK(x14.to_string() == "SE(2)");
    CHECK(x15.to_string() == "SO(2)");
    CHECK(x16.to_string() == "SP(2, R)");
    CHECK(x17.to_string() == "SU(2)");
    CHECK(x1m8.to_string() == "C^2");
    CHECK(x1m7.to_string() == "GL(2, C)");
    CHECK(x1m6.to_string() == "GL(2, R)");
    CHECK(x1m5.to_string() == "R^2");
    CHECK(x1m4.to_string() == "SE(2)");
    CHECK(x1m3.to_string() == "SO(2)");
    CHECK(x1m2.to_string() == "SP(2, R)");
    CHECK(x1m1.to_string() == "SU(2)");

    // Out of bounds
    const GLC x18 = std::get<GLC>(x1[8]);
    const GLC x1m9 = std::get<GLC>(x1[-9]);

    CHECK(x18.to_string() == "GL(0, C)");
    CHECK(x1m9.to_string() == "GL(0, C)");
}
