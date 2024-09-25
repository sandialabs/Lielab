#include <cmath>
#include <Lielab.hpp>

constexpr double PI = Lielab::constants::PI<double>;
const Lielab::domain::SO DCMId(3);
const Lielab::domain::SO DCMrotx(Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(0,3)));
const Lielab::domain::SO DCMroty(Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(1,3)));
const Lielab::domain::SO DCMrotz(Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(2,3)));
const double some_angle = PI/2.0*5.0/7.0;

TEST_CASE("eanglebody123_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody123_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody123_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody123_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody123_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody123_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody123_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody231_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody231_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody231_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody231_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody231_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody231_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody231_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody312_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody312_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody312_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody312_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody312_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody312_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody312_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody132_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody132_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody132_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody132_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody132_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody132_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody132_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody213_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody213_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody213_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody213_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody213_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody213_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody213_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody321_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody321_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody321_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody321_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody321_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody321_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody321_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody121_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody121_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody121_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody121_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody121_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody121_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody121_to_dcm<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody121_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody131_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody131_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody131_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody131_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody131_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody131_to_dcm<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody131_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody131_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody212_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody212_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody212_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody212_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody212_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody212_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody212_to_dcm<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody212_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody232_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody232_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody232_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody232_to_dcm<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody232_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody232_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody232_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody232_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody313_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody313_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody313_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody313_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody313_to_dcm<double>(-PI/2.0, -PI/2.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody313_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody313_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody313_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody323_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody323_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglebody323_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglebody323_to_dcm<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglebody323_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody323_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglebody323_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglebody323_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace123_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace123_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace123_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace123_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace123_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace123_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace123_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace231_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace231_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace231_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace231_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace231_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace231_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace231_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace312_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace231_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace312_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace312_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace312_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace312_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace312_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace132_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace132_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace132_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace132_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace132_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace132_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace132_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace213_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace213_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace213_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace213_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace213_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace213_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace213_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace321_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace321_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace321_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace321_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace321_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace321_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace321_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace121_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace121_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace121_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace121_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace121_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace121_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace121_to_dcm<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace121_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace131_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace131_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace131_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace131_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace131_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace131_to_dcm<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace131_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace131_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace212_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace212_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace212_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace212_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace212_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace212_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace212_to_dcm<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace212_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace232_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace232_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace232_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace232_to_dcm<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace232_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace232_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace232_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace232_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace313_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace313_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace313_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace313_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace313_to_dcm<double>(-PI/2.0, PI/2.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace313_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace313_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace313_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace323_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace323_to_dcm function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::transform::eanglespace323_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::transform::eanglespace323_to_dcm<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::transform::eanglespace323_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace323_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::transform::eanglespace323_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::transform::eanglespace323_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}
