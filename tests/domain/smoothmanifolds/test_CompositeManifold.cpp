#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

#include <iostream>

Lielab::domain::CompositeManifold _make_cmanifold()
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

    const Lielab::domain::cn ycn1 = Lielab::domain::cn::from_vector({44.0, 45.0, 46.0, 47.0});
    const Lielab::domain::glc yglc1 = Lielab::domain::glc::from_vector({48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0});
    const Lielab::domain::glr yglr1 = Lielab::domain::glr::from_vector({56.0, 57.0, 58.0, 59.0});
    const Lielab::domain::rn yrn1 = Lielab::domain::rn::from_vector({60.0, 61.0});
    const Lielab::domain::se yse1 = Lielab::domain::se::from_vector({62.0, 63.0, 64.0});
    const Lielab::domain::so yso1 = Lielab::domain::so::from_vector({65.0, 66.0, 67.0});
    const Lielab::domain::sp ysp1 = Lielab::domain::sp::from_vector({68.0, 69.0, 70.0});
    const Lielab::domain::su ysu1 = Lielab::domain::su::from_vector({71.0, 72.0, 73.0});

    return CompositeManifold({yCN1, yGLC1, yGLR1, yRN1, ySE1, ySO1, ySP1, ySU1, ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});
}

TEST_CASE("CompositeManifold to_string", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeManifold y1 = _make_cmanifold();

    CHECK(y1.to_string() == "C^2 x GL(2, C) x GL(2, R) x R^2 x SE(2) x SO(2) x SP(2, R) x SU(2) x c^2 x gl(2, C) x gl(2, R) x r^2 x se(2) x so(3) x sp(2, R) x su(2)");
}

TEST_CASE("CompositeManifold main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeManifold xblank = CompositeManifold();
    CHECK(xblank.get_dimension() == 0);

    const CompositeManifold x0 = CompositeManifold(0);
    CHECK(x0.get_dimension() == 0);
    const CompositeManifold x1 = CompositeManifold(1);
    CHECK(x1.get_dimension() == 2);
    const CompositeManifold x10 = CompositeManifold(10);
    CHECK(x10.get_dimension() == 200);
}

TEST_CASE("CompositeManifold from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeManifold x0 = CompositeManifold::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXcd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const CompositeManifold x1 = CompositeManifold::from_shape(1);
    CHECK(x1.get_dimension() == 2);
    const Eigen::MatrixXcd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const CompositeManifold x2 = CompositeManifold::from_shape(2);
    CHECK(x2.get_dimension() == 8);
    const Eigen::MatrixXcd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("CompositeManifold get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    CompositeManifold zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

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

    CompositeManifold y1 = _make_cmanifold();
    
    const std::vector<size_t> dims = y1.get_dimensions();
    REQUIRE(dims.size() == 16);
    CHECK(dims[0] == 4);
    CHECK(dims[1] == 8);
    CHECK(dims[2] == 4);
    CHECK(dims[3] == 2);
    CHECK(dims[4] == 3);
    CHECK(dims[5] == 1);
    CHECK(dims[6] == 3);
    CHECK(dims[7] == 3);
    CHECK(dims[8] == 4);
    CHECK(dims[9] == 8);
    CHECK(dims[10] == 4);
    CHECK(dims[11] == 2);
    CHECK(dims[12] == 3);
    CHECK(dims[13] == 3);
    CHECK(dims[14] == 3);
    CHECK(dims[15] == 3);

    CHECK(y1.get_dimension() == 58);
}

TEST_CASE("CompositeManifold get_size", "[domain]")
{
    using namespace Lielab::domain;

    CompositeManifold zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

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

    CompositeManifold y1 = _make_cmanifold();
    
    const std::vector<size_t> sizes = y1.get_sizes();
    REQUIRE(sizes.size() == 16);
    CHECK(sizes[0] == 4);
    CHECK(sizes[1] == 8);
    CHECK(sizes[2] == 4);
    CHECK(sizes[3] == 2);
    CHECK(sizes[4] == 9);
    CHECK(sizes[5] == 4);
    CHECK(sizes[6] == 4);
    CHECK(sizes[7] == 8);
    CHECK(sizes[8] == 4);
    CHECK(sizes[9] == 8);
    CHECK(sizes[10] == 4);
    CHECK(sizes[11] == 2);
    CHECK(sizes[12] == 3);
    CHECK(sizes[13] == 3);
    CHECK(sizes[14] == 3);
    CHECK(sizes[15] == 3);

    CHECK(y1.get_size() == 73);
}

TEST_CASE("CompositeManifold serialize/unserialize", "[domain]")
{
    /*!
    * Tests the serialize/unserialize operation.
    */

    using namespace Lielab::domain;

    CompositeManifold y1 = _make_cmanifold();
    y1.unserialize({73.0, 72.0, 71.0, 70.0, 69.0, 68.0, 67.0, 66.0, 65.0, 64.0, 63.0, 62.0, 61.0, 60.0,
                    59.0, 58.0, 57.0, 56.0, 55.0, 54.0, 53.0, 52.0, 51.0, 50.0, 49.0, 48.0, 47.0, 46.0,
                    45.0, 44.0, 43.0, 42.0, 41.0, 40.0, 39.0, 38.0, 37.0, 36.0, 35.0, 34.0, 33.0, 32.0,
                    31.0, 30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0, 20.0, 19.0, 18.0,
                    17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0,
                    2.0, 1.0});
    
    const Eigen::VectorXd y1bar = y1.serialize();
    REQUIRE(y1bar.size() == 73);
    CHECK(y1bar(0) == 73.0);
    CHECK(y1bar(1) == 72.0);
    CHECK(y1bar(2) == 71.0);
    CHECK(y1bar(3) == 70.0);
    CHECK(y1bar(4) == 69.0);
    CHECK(y1bar(5) == 68.0);
    CHECK(y1bar(6) == 67.0);
    CHECK(y1bar(7) == 66.0);
    CHECK(y1bar(8) == 65.0);
    CHECK(y1bar(9) == 64.0);
    CHECK(y1bar(10) == 63.0);
    CHECK(y1bar(11) == 62.0);
    CHECK(y1bar(12) == 61.0);
    CHECK(y1bar(13) == 60.0);
    CHECK(y1bar(14) == 59.0);
    CHECK(y1bar(15) == 58.0);
    CHECK(y1bar(16) == 57.0);
    CHECK(y1bar(17) == 56.0);
    CHECK(y1bar(18) == 55.0);
    CHECK(y1bar(19) == 54.0);
    CHECK(y1bar(20) == 53.0);
    CHECK(y1bar(21) == 52.0);
    CHECK(y1bar(22) == 51.0);
    CHECK(y1bar(23) == 50.0);
    CHECK(y1bar(24) == 49.0);
    CHECK(y1bar(25) == 48.0);
    CHECK(y1bar(26) == 47.0);
    CHECK(y1bar(27) == 46.0);
    CHECK(y1bar(28) == 45.0);
    CHECK(y1bar(29) == 44.0);
    CHECK(y1bar(30) == 43.0);
    CHECK(y1bar(31) == 42.0);
    CHECK(y1bar(32) == 41.0);
    CHECK(y1bar(33) == 40.0);
    CHECK(y1bar(34) == 39.0);
    CHECK(y1bar(35) == 38.0);
    CHECK(y1bar(36) == 37.0);
    CHECK(y1bar(37) == 36.0);
    CHECK(y1bar(38) == 35.0);
    CHECK(y1bar(39) == 34.0);
    CHECK(y1bar(40) == 33.0);
    CHECK(y1bar(41) == 32.0);
    CHECK(y1bar(42) == 31.0);
    CHECK(y1bar(43) == 30.0);
    CHECK(y1bar(44) == 29.0);
    CHECK(y1bar(45) == 28.0);
    CHECK(y1bar(46) == 27.0);
    CHECK(y1bar(47) == 26.0);
    CHECK(y1bar(48) == 25.0);
    CHECK(y1bar(49) == 24.0);
    CHECK(y1bar(50) == 23.0);
    CHECK(y1bar(51) == 22.0);
    CHECK(y1bar(52) == 21.0);
    CHECK(y1bar(53) == 20.0);
    CHECK(y1bar(54) == 19.0);
    CHECK(y1bar(55) == 18.0);
    CHECK(y1bar(56) == 17.0);
    CHECK(y1bar(57) == 16.0);
    CHECK(y1bar(58) == 15.0);
    CHECK(y1bar(59) == 14.0);
    CHECK(y1bar(60) == 13.0);
    CHECK(y1bar(61) == 12.0);
    CHECK(y1bar(62) == 11.0);
    CHECK(y1bar(63) == 10.0);
    CHECK(y1bar(64) == 9.0);
    CHECK(y1bar(65) == 8.0);
    CHECK(y1bar(66) == 7.0);
    CHECK(y1bar(67) == 6.0);
    CHECK(y1bar(68) == 5.0);
    CHECK(y1bar(69) == 4.0);
    CHECK(y1bar(70) == 3.0);
    CHECK(y1bar(71) == 2.0);
    CHECK(y1bar(72) == 1.0);
}

TEST_CASE("CompositeManifold get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    const CompositeManifold y1 = _make_cmanifold();

    const Eigen::MatrixXcd y1hat = y1.get_matrix();

    REQUIRE(y1hat.rows() == 39);
    REQUIRE(y1hat.cols() == 39);
    
    // TODO: Check the 0's
    // TODO: Check each submatrix

}

TEST_CASE("CompositeManifold operator()", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeManifold y1 = _make_cmanifold();

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
    CHECK(y1(-39, -39) == std::complex<double>(1.0, 0.0));
    CHECK(y1(-39, -38) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-39, -37) == std::complex<double>(1.0, 2.0));
    CHECK(y1(-38, -39) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-38, -38) == std::complex<double>(1.0, 0.0));
    CHECK(y1(-38, -37) == std::complex<double>(3.0, 4.0));
    CHECK(y1(-37, -39) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-37, -38) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-37, -37) == std::complex<double>(1.0, 0.0));

    // In bounds GLC component
    CHECK(y1(3, 3) == std::complex<double>(5.0, 6.0));
    CHECK(y1(3, 4) == std::complex<double>(7.0, 8.0));
    CHECK(y1(4, 3) == std::complex<double>(9.0, 10.0));
    CHECK(y1(4, 4) == std::complex<double>(11.0, 12.0));
    CHECK(y1(-36, -36) == std::complex<double>(5.0, 6.0));
    CHECK(y1(-36, -35) == std::complex<double>(7.0, 8.0));
    CHECK(y1(-35, -36) == std::complex<double>(9.0, 10.0));
    CHECK(y1(-35, -35) == std::complex<double>(11.0, 12.0));

    // In bounds GLR component
    CHECK(y1(5, 5) == std::complex<double>(13.0, 0.0));
    CHECK(y1(5, 6) == std::complex<double>(14.0, 0.0));
    CHECK(y1(6, 5) == std::complex<double>(15.0, 0.0));
    CHECK(y1(6, 6) == std::complex<double>(16.0, 0.0));
    CHECK(y1(-34, -34) == std::complex<double>(13.0, 0.0));
    CHECK(y1(-34, -33) == std::complex<double>(14.0, 0.0));
    CHECK(y1(-33, -34) == std::complex<double>(15.0, 0.0));
    CHECK(y1(-33, -33) == std::complex<double>(16.0, 0.0));

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
    CHECK(y1(-32, -32) == std::complex<double>(1.0, 0.0));
    CHECK(y1(-32, -31) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-32, -30) == std::complex<double>(17.0, 0.0));
    CHECK(y1(-31, -32) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-31, -31) == std::complex<double>(1.0, 0.0));
    CHECK(y1(-31, -30) == std::complex<double>(18.0, 0.0));
    CHECK(y1(-30, -32) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-30, -31) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-30, -30) == std::complex<double>(1.0, 0.0));

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
    CHECK(y1(-29, -29) == std::complex<double>(19.0, 0.0));
    CHECK(y1(-29, -28) == std::complex<double>(20.0, 0.0));
    CHECK(y1(-29, -27) == std::complex<double>(21.0, 0.0));
    CHECK(y1(-28, -29) == std::complex<double>(22.0, 0.0));
    CHECK(y1(-28, -28) == std::complex<double>(23.0, 0.0));
    CHECK(y1(-28, -27) == std::complex<double>(24.0, 0.0));
    CHECK(y1(-27, -29) == std::complex<double>(25.0, 0.0));
    CHECK(y1(-27, -28) == std::complex<double>(26.0, 0.0));
    CHECK(y1(-27, -27) == std::complex<double>(27.0, 0.0));

    // In bounds SO component
    CHECK(y1(13, 13) == std::complex<double>(28.0, 0.0));
    CHECK(y1(13, 14) == std::complex<double>(29.0, 0.0));
    CHECK(y1(14, 13) == std::complex<double>(30.0, 0.0));
    CHECK(y1(14, 14) == std::complex<double>(31.0, 0.0));
    CHECK(y1(-26, -26) == std::complex<double>(28.0, 0.0));
    CHECK(y1(-26, -25) == std::complex<double>(29.0, 0.0));
    CHECK(y1(-25, -26) == std::complex<double>(30.0, 0.0));
    CHECK(y1(-25, -25) == std::complex<double>(31.0, 0.0));

    // In bounds SP component
    CHECK(y1(15, 15) == std::complex<double>(32.0, 0.0));
    CHECK(y1(15, 16) == std::complex<double>(33.0, 0.0));
    CHECK(y1(16, 15) == std::complex<double>(34.0, 0.0));
    CHECK(y1(16, 16) == std::complex<double>(35.0, 0.0));
    CHECK(y1(-24, -24) == std::complex<double>(32.0, 0.0));
    CHECK(y1(-24, -23) == std::complex<double>(33.0, 0.0));
    CHECK(y1(-23, -24) == std::complex<double>(34.0, 0.0));
    CHECK(y1(-23, -23) == std::complex<double>(35.0, 0.0));

    // In bounds SU component
    CHECK(y1(17, 17) == std::complex<double>(36.0, 37.0));
    CHECK(y1(17, 18) == std::complex<double>(38.0, 39.0));
    CHECK(y1(18, 17) == std::complex<double>(40.0, 41.0));
    CHECK(y1(18, 18) == std::complex<double>(42.0, 43.0));
    CHECK(y1(-22, -22) == std::complex<double>(36.0, 37.0));
    CHECK(y1(-22, -21) == std::complex<double>(38.0, 39.0));
    CHECK(y1(-21, -22) == std::complex<double>(40.0, 41.0));
    CHECK(y1(-21, -21) == std::complex<double>(42.0, 43.0));

    // In bounds cn component
    CHECK(y1(19, 19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(19, 20) == std::complex<double>(0.0, 0.0));
    CHECK(y1(19, 21) == std::complex<double>(44.0, 45.0));
    CHECK(y1(20, 19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(20, 20) == std::complex<double>(0.0, 0.0));
    CHECK(y1(20, 21) == std::complex<double>(46.0, 47.0));
    CHECK(y1(21, 19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(21, 20) == std::complex<double>(0.0, 0.0));
    CHECK(y1(21, 21) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-20, -20) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-20, -19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-20, -18) == std::complex<double>(44.0, 45.0));
    CHECK(y1(-19, -20) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-19, -19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-19, -18) == std::complex<double>(46.0, 47.0));
    CHECK(y1(-18, -20) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-18, -19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-18, -18) == std::complex<double>(0.0, 0.0));

    // In bounds glc component
    CHECK(y1(22, 22) == std::complex<double>(48.0, 49.0));
    CHECK(y1(22, 23) == std::complex<double>(50.0, 51.0));
    CHECK(y1(23, 22) == std::complex<double>(52.0, 53.0));
    CHECK(y1(23, 23) == std::complex<double>(54.0, 55.0));
    CHECK(y1(-17, -17) == std::complex<double>(48.0, 49.0));
    CHECK(y1(-17, -16) == std::complex<double>(50.0, 51.0));
    CHECK(y1(-16, -17) == std::complex<double>(52.0, 53.0));
    CHECK(y1(-16, -16) == std::complex<double>(54.0, 55.0));

    // In bounds glr component
    CHECK(y1(24, 24) == std::complex<double>(56.0, 0.0));
    CHECK(y1(24, 25) == std::complex<double>(57.0, 0.0));
    CHECK(y1(25, 24) == std::complex<double>(58.0, 0.0));
    CHECK(y1(25, 25) == std::complex<double>(59.0, 0.0));
    CHECK(y1(-15, -15) == std::complex<double>(56.0, 0.0));
    CHECK(y1(-15, -14) == std::complex<double>(57.0, 0.0));
    CHECK(y1(-14, -15) == std::complex<double>(58.0, 0.0));
    CHECK(y1(-14, -14) == std::complex<double>(59.0, 0.0));

    // In bounds rn component
    CHECK(y1(26, 26) == std::complex<double>(0.0, 0.0));
    CHECK(y1(26, 27) == std::complex<double>(0.0, 0.0));
    CHECK(y1(26, 28) == std::complex<double>(60.0, 0.0));
    CHECK(y1(27, 26) == std::complex<double>(0.0, 0.0));
    CHECK(y1(27, 27) == std::complex<double>(0.0, 0.0));
    CHECK(y1(27, 28) == std::complex<double>(61.0, 0.0));
    CHECK(y1(28, 26) == std::complex<double>(0.0, 0.0));
    CHECK(y1(28, 27) == std::complex<double>(0.0, 0.0));
    CHECK(y1(28, 28) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-13, -13) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-13, -12) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-13, -11) == std::complex<double>(60.0, 0.0));
    CHECK(y1(-12, -13) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-12, -12) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-12, -11) == std::complex<double>(61.0, 0.0));
    CHECK(y1(-11, -13) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-11, -12) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-11, -11) == std::complex<double>(0.0, 0.0));

    // In bounds se component
    CHECK(y1(29, 29) == std::complex<double>(0.0, 0.0));
    CHECK(y1(29, 30) == std::complex<double>(-64.0, 0.0));
    CHECK(y1(29, 31) == std::complex<double>(62.0, 0.0));
    CHECK(y1(30, 29) == std::complex<double>(64.0, 0.0));
    CHECK(y1(30, 30) == std::complex<double>(0.0, 0.0));
    CHECK(y1(30, 31) == std::complex<double>(63.0, 0.0));
    CHECK(y1(31, 29) == std::complex<double>(0.0, 0.0));
    CHECK(y1(31, 30) == std::complex<double>(0.0, 0.0));
    CHECK(y1(31, 31) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-10, -10) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-10, -9) == std::complex<double>(-64.0, 0.0));
    CHECK(y1(-10, -8) == std::complex<double>(62.0, 0.0));
    CHECK(y1(-9, -10) == std::complex<double>(64.0, 0.0));
    CHECK(y1(-9, -9) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-9, -8) == std::complex<double>(63.0, 0.0));
    CHECK(y1(-8, -10) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-8, -9) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-8, -8) == std::complex<double>(0.0, 0.0));

    // In bounds so component
    CHECK(y1(32, 32) == std::complex<double>(0.0, 0.0));
    CHECK(y1(32, 33) == std::complex<double>(-67.0, 0.0));
    CHECK(y1(32, 34) == std::complex<double>(66.0, 0.0));
    CHECK(y1(33, 32) == std::complex<double>(67.0, 0.0));
    CHECK(y1(33, 33) == std::complex<double>(0.0, 0.0));
    CHECK(y1(33, 34) == std::complex<double>(-65.0, 0.0));
    CHECK(y1(34, 32) == std::complex<double>(-66.0, 0.0));
    CHECK(y1(34, 33) == std::complex<double>(65.0, 0.0));
    CHECK(y1(34, 34) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-7, -7) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-7, -6) == std::complex<double>(-67.0, 0.0));
    CHECK(y1(-7, -5) == std::complex<double>(66.0, 0.0));
    CHECK(y1(-6, -7) == std::complex<double>(67.0, 0.0));
    CHECK(y1(-6, -6) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-6, -5) == std::complex<double>(-65.0, 0.0));
    CHECK(y1(-5, -7) == std::complex<double>(-66.0, 0.0));
    CHECK(y1(-5, -6) == std::complex<double>(65.0, 0.0));
    CHECK(y1(-5, -5) == std::complex<double>(0.0, 0.0));

    // In bounds sp component
    CHECK(y1(35, 35) == std::complex<double>(68.0, 0.0));
    CHECK(y1(35, 36) == std::complex<double>(69.0, 0.0));
    CHECK(y1(36, 35) == std::complex<double>(70.0, 0.0));
    CHECK(y1(36, 36) == std::complex<double>(-68.0, 0.0));
    CHECK(y1(-4, -4) == std::complex<double>(68.0, 0.0));
    CHECK(y1(-4, -3) == std::complex<double>(69.0, 0.0));
    CHECK(y1(-3, -4) == std::complex<double>(70.0, 0.0));
    CHECK(y1(-3, -3) == std::complex<double>(-68.0, 0.0));

    // In bounds su component
    CHECK(y1(37, 37) == std::complex<double>(0.0, 73.0));
    CHECK(y1(37, 38) == std::complex<double>(-72.0, 71.0));
    CHECK(y1(38, 37) == std::complex<double>(72.0, 71.0));
    CHECK(y1(38, 38) == std::complex<double>(0.0, -73.0));
    CHECK(y1(-2, -2) == std::complex<double>(0.0, 73.0));
    CHECK(y1(-2, -1) == std::complex<double>(-72.0, 71.0));
    CHECK(y1(-1, -2) == std::complex<double>(72.0, 71.0));
    CHECK(y1(-1, -1) == std::complex<double>(0.0, -73.0));

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
    CHECK(std::isnan(y1(0, -40).real()));
    CHECK(std::isnan(y1(0, -40).imag()));
    CHECK(std::isnan(y1(-40, 0).real()));
    CHECK(std::isnan(y1(-40, 0).imag()));
    CHECK(std::isnan(y1(-40, -40).real()));
    CHECK(std::isnan(y1(-40, -40).imag()));
    CHECK(std::isnan(y1(0, 39).real()));
    CHECK(std::isnan(y1(0, 39).imag()));
    CHECK(std::isnan(y1(39, 0).real()));
    CHECK(std::isnan(y1(39, 0).imag()));
    CHECK(std::isnan(y1(39, 39).real()));
    CHECK(std::isnan(y1(39, 39).imag()));
}

TEST_CASE("CompositeManifold operator[]", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeManifold x1 = _make_cmanifold();
    
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
    const cn x18 = std::get<cn>(x1[8]);
    const glc x19 = std::get<glc>(x1[9]);
    const glr x110 = std::get<glr>(x1[10]);
    const rn x111 = std::get<rn>(x1[11]);
    const se x112 = std::get<se>(x1[12]);
    const so x113 = std::get<so>(x1[13]);
    const sp x114 = std::get<sp>(x1[14]);
    const su x115 = std::get<su>(x1[15]);
    const CN x1m16 = std::get<CN>(x1[-16]);
    const GLC x1m15 = std::get<GLC>(x1[-15]);
    const GLR x1m14 = std::get<GLR>(x1[-14]);
    const RN x1m13 = std::get<RN>(x1[-13]);
    const SE x1m12 = std::get<SE>(x1[-12]);
    const SO x1m11 = std::get<SO>(x1[-11]);
    const SP x1m10 = std::get<SP>(x1[-10]);
    const SU x1m9 = std::get<SU>(x1[-9]);
    const cn x1m8 = std::get<cn>(x1[-8]);
    const glc x1m7 = std::get<glc>(x1[-7]);
    const glr x1m6 = std::get<glr>(x1[-6]);
    const rn x1m5 = std::get<rn>(x1[-5]);
    const se x1m4 = std::get<se>(x1[-4]);
    const so x1m3 = std::get<so>(x1[-3]);
    const sp x1m2 = std::get<sp>(x1[-2]);
    const su x1m1 = std::get<su>(x1[-1]);
    

    CHECK(x10.to_string() == "C^2");
    CHECK(x11.to_string() == "GL(2, C)");
    CHECK(x12.to_string() == "GL(2, R)");
    CHECK(x13.to_string() == "R^2");
    CHECK(x14.to_string() == "SE(2)");
    CHECK(x15.to_string() == "SO(2)");
    CHECK(x16.to_string() == "SP(2, R)");
    CHECK(x17.to_string() == "SU(2)");
    CHECK(x18.to_string() == "c^2");
    CHECK(x19.to_string() == "gl(2, C)");
    CHECK(x110.to_string() == "gl(2, R)");
    CHECK(x111.to_string() == "r^2");
    CHECK(x112.to_string() == "se(2)");
    CHECK(x113.to_string() == "so(3)");
    CHECK(x114.to_string() == "sp(2, R)");
    CHECK(x115.to_string() == "su(2)");
    CHECK(x1m16.to_string() == "C^2");
    CHECK(x1m15.to_string() == "GL(2, C)");
    CHECK(x1m14.to_string() == "GL(2, R)");
    CHECK(x1m13.to_string() == "R^2");
    CHECK(x1m12.to_string() == "SE(2)");
    CHECK(x1m11.to_string() == "SO(2)");
    CHECK(x1m10.to_string() == "SP(2, R)");
    CHECK(x1m9.to_string() == "SU(2)");
    CHECK(x1m8.to_string() == "c^2");
    CHECK(x1m7.to_string() == "gl(2, C)");
    CHECK(x1m6.to_string() == "gl(2, R)");
    CHECK(x1m5.to_string() == "r^2");
    CHECK(x1m4.to_string() == "se(2)");
    CHECK(x1m3.to_string() == "so(3)");
    CHECK(x1m2.to_string() == "sp(2, R)");
    CHECK(x1m1.to_string() == "su(2)");

    // Out of bounds
    const glc x116 = std::get<glc>(x1[16]);
    const glc x1m17 = std::get<glc>(x1[-17]);

    CHECK(x116.to_string() == "gl(0, C)");
    CHECK(x1m17.to_string() == "gl(0, C)");
}
