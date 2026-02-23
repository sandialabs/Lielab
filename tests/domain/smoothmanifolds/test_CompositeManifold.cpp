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
    ySE1.unserialize({19.0, 20.0, 21.0, 22.0, 23.0, 24.0});
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
    CHECK(x0.point.size() == 0);
    const CompositeManifold x1 = CompositeManifold(1);
    CHECK(x1.point.size() == 1);
    const CompositeManifold x10 = CompositeManifold(10);
    CHECK(x10.point.size() == 10);
}

TEST_CASE("CompositeManifold get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    CompositeManifold y1 = _make_cmanifold();
    
    const std::vector<int> dims = y1.get_dimensions();
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

    CompositeManifold y1 = _make_cmanifold();
    
    const std::vector<int> sizes = y1.get_sizes();
    REQUIRE(sizes.size() == 16);
    CHECK(sizes[0] == 4);
    CHECK(sizes[1] == 8);
    CHECK(sizes[2] == 4);
    CHECK(sizes[3] == 2);
    CHECK(sizes[4] == 6);
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

    CHECK(y1.get_size() == 70);
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
                    17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0});
    
    const Eigen::VectorXd y1bar = y1.serialize();
    REQUIRE(y1bar.size() == 70);
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
}

TEST_CASE("CompositeManifold operator[]", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeManifold x1 = _make_cmanifold();
    
    const CN x10 = x1[0];
    const GLC x11 = x1[1];
    const GLR x12 = x1[2];
    const RN x13 = x1[3];
    const SE x14 = x1[4];
    const SO x15 = x1[5];
    const SP x16 = x1[6];
    const SU x17 = x1[7];
    const cn x18 = x1[8];
    const glc x19 = x1[9];
    const glr x110 = x1[10];
    const rn x111 = x1[11];
    const se x112 = x1[12];
    const so x113 = x1[13];
    const sp x114 = x1[14];
    const su x115 = x1[15];
    const CN x1m16 = x1[-16];
    const GLC x1m15 = x1[-15];
    const GLR x1m14 = x1[-14];
    const RN x1m13 = x1[-13];
    const SE x1m12 = x1[-12];
    const SO x1m11 = x1[-11];
    const SP x1m10 = x1[-10];
    const SU x1m9 = x1[-9];
    const cn x1m8 = x1[-8];
    const glc x1m7 = x1[-7];
    const glr x1m6 = x1[-6];
    const rn x1m5 = x1[-5];
    const se x1m4 = x1[-4];
    const so x1m3 = x1[-3];
    const sp x1m2 = x1[-2];
    const su x1m1 = x1[-1];
    

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
    CHECK_THROWS(x1[16]);
    CHECK_THROWS(x1[-17]);
}
