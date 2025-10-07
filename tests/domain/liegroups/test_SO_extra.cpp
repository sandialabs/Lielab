#include <cmath>
#include <numbers>
#include <Lielab.hpp>
#include <catch2/catch_all.hpp>

#include "../../test_utils.hpp"

constexpr double PI = std::numbers::pi_v<double>;
const Lielab::domain::SO DCMId(3);
const Lielab::domain::SO DCMrotx(Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(0,3)));
const Lielab::domain::SO DCMroty(Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(1,3)));
const Lielab::domain::SO DCMrotz(Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(2,3)));
const double some_angle = PI/2.0*5.0/7.0;

TEST_CASE("SO", "[domain]")
{
    /*!
    * Tests SO against well-known identities.
    */

    for (size_t shape = 2; shape <= 4; shape++)
    {
        const size_t D = Lielab::domain::so::basis(0, shape).get_dimension();

        // Construct the SO elements
        std::vector<Lielab::domain::SO> elements;
        for (size_t ii = 0; ii < D; ii++)
        {
            elements.push_back(Lielab::functions::exp(Lielab::domain::so::basis(ii, shape)));
        }

        const Lielab::domain::SO identity(shape);

        is_group<Lielab::domain::SO>(elements, identity);
    }
}

TEST_CASE("from_SU2", "[domain]")
{
    /*!
    * Tests the from_SU2 function
    */

    // Build 90 degree rotations in x, y, and z
    // TODO: simplify with basis()
    Lielab::domain::su u(2), v(2), w(2);
    Eigen::Vector3d xx, yy, zz;
    xx << 0.0, 0.0, 0.5;
    yy << 0.0, 0.5, 0.0;
    zz << 0.5, 0.0, 0.0;
    u.set_vector(xx);
    v.set_vector(yy);
    w.set_vector(zz);
    Lielab::domain::SU qx = Lielab::functions::exp(PI/2.0*u);
    Lielab::domain::SU qy = Lielab::functions::exp(PI/2.0*v);
    Lielab::domain::SU qz = Lielab::functions::exp(PI/2.0*w);

    // Test 90 degree x rotation
    Lielab::domain::SO rx = Lielab::domain::SO::from_SU2(qx);

    assert_domain(rx, DCMrotx);

    // Test 90 degree y rotation
    Lielab::domain::SO ry = Lielab::domain::SO::from_SU2(qy);

    assert_domain(ry, DCMroty);

    // Test 90 degree z rotation
    Lielab::domain::SO rz = Lielab::domain::SO::from_SU2(qz);

    assert_domain(rz, DCMrotz);
}

TEST_CASE("from_eulerangles_body123", "[domain]")
{
    /*!
    * Tests the eanglebody123 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body123<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body123<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body123<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body123<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body123<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body231", "[domain]")
{
    /*!
    * Tests the eanglebody231 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body231<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body231<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body231<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body231<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body231<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body312", "[domain]")
{
    /*!
    * Tests the eanglebody312 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body312<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body312<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body312<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body312<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body312<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body132", "[domain]")
{
    /*!
    * Tests the eanglebody132 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body132<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body132<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body132<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body132<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body132<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body213", "[domain]")
{
    /*!
    * Tests the eanglebody213 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body213<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body213<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body213<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body213<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body213<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body321", "[domain]")
{
    /*!
    * Tests the eanglebody321 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body321<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body321<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body321<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body321<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body321<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body121", "[domain]")
{
    /*!
    * Tests the eanglebody121 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body121<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body121<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body121<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body121<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body121<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body121<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body131", "[domain]")
{
    /*!
    * Tests the eanglebody131 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body131<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body131<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body131<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body131<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body131<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body131<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body212", "[domain]")
{
    /*!
    * Tests the eanglebody212 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body212<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body212<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body212<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body212<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body212<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body212<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body232", "[domain]")
{
    /*!
    * Tests the eanglebody232 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body232<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body232<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body232<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body232<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body232<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body232<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body313", "[domain]")
{
    /*!
    * Tests the eanglebody313 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body313<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body313<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body313<double>(-PI/2.0, -PI/2.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body313<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body313<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body313<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_body323", "[domain]")
{
    /*!
    * Tests the eanglebody323 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_body323<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body323<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body323<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body323<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body323<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body323<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space123", "[domain]")
{
    /*!
    * Tests the eanglespace123 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space123<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space123<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space123<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space123<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space123<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space231", "[domain]")
{
    /*!
    * Tests the eanglespace231 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space231<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space231<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space231<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space231<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space231<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space312", "[domain]")
{
    /*!
    * Tests the eanglespace231 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space312<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space312<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space312<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space312<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space312<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space132", "[domain]")
{
    /*!
    * Tests the eanglespace132 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space132<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space132<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space132<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space132<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space132<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space213", "[domain]")
{
    /*!
    * Tests the eanglespace213 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space213<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space213<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space213<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space213<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space213<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space321", "[domain]")
{
    /*!
    * Tests the eanglespace321 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space321<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space321<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space321<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space321<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space321<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space121", "[domain]")
{
    /*!
    * Tests the eanglespace121 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space121<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space121<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space121<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space121<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space121<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space121<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space131", "[domain]")
{
    /*!
    * Tests the eanglespace131 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space131<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space131<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space131<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space131<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space131<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space131<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space212", "[domain]")
{
    /*!
    * Tests the eanglespace212 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space212<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space212<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space212<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space212<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space212<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space212<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space232", "[domain]")
{
    /*!
    * Tests the eanglespace232 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space232<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space232<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space232<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space232<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space232<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space232<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space313", "[domain]")
{
    /*!
    * Tests the eanglespace313 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space313<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space313<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space313<double>(-PI/2.0, PI/2.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space313<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space313<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space313<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("from_eulerangles_space323", "[domain]")
{
    /*!
    * Tests the eanglespace323 function.
    */

    // Identity
    Lielab::domain::SO dcm = Lielab::domain::SO::from_eulerangles_space323<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space323<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space323<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space323<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space323<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space323<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_matrix().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_matrix()*dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("to_eulerangles_body123", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body123 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body123<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body123<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body123<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body123<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body123<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body123<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body123(0.5, PI/2.0, 0.5).to_eulerangles_body123<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body123(1.0, PI/2.0, 0.0).to_eulerangles_body123<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body123(1.0, -PI/2.0, 0.0).to_eulerangles_body123<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body123(some_angle, some_angle, some_angle).to_eulerangles_body123<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body231", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body231 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body231<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body231<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body231<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body231<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body231<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body231<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body231(0.5, PI/2.0, 0.5).to_eulerangles_body231<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body231(1.0, PI/2.0, 0.0).to_eulerangles_body231<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body231(1.0, -PI/2.0, 0.0).to_eulerangles_body231<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body231(some_angle, some_angle, some_angle).to_eulerangles_body231<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body312", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body312 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body312<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body312<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body312<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body312<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body312<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body312<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body312(0.5, PI/2.0, 0.5).to_eulerangles_body312<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body312(1.0, PI/2.0, 0.0).to_eulerangles_body312<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body312(1.0, -PI/2.0, 0.0).to_eulerangles_body312<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body312(some_angle, some_angle, some_angle).to_eulerangles_body312<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body132", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body132 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body132<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body132<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body132<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body132<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body132<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body132<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body132(0.5, PI/2.0, 0.5).to_eulerangles_body132<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body132(1.0, PI/2.0, 0.0).to_eulerangles_body132<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body132(1.0, -PI/2.0, 0.0).to_eulerangles_body132<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body132(some_angle, some_angle, some_angle).to_eulerangles_body132<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body213", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body213 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body213<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body213<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body213<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body213<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body213<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body213<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body213(0.5, PI/2.0, 0.5).to_eulerangles_body213<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body213(1.0, PI/2.0, 0.0).to_eulerangles_body213<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body213(1.0, -PI/2.0, 0.0).to_eulerangles_body213<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body213(some_angle, some_angle, some_angle).to_eulerangles_body213<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body321", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body321 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body321<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body321<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body321<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body321<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body321<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body321<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body321(0.5, PI/2.0, 0.5).to_eulerangles_body321<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body321(1.0, PI/2.0, 0.0).to_eulerangles_body321<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body321(1.0, -PI/2.0, 0.0).to_eulerangles_body321<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body321(some_angle, some_angle, some_angle).to_eulerangles_body321<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body121", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body121 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body121<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body121<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body121<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body121<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 1));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body121<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body121<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(-PI/2.0, 1));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body121(0.5, 0.0, 0.5).to_eulerangles_body121<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body121(1.0, 0.0, 0.0).to_eulerangles_body121<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body121(1.0, PI, 0.0).to_eulerangles_body121<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body121(some_angle, some_angle, some_angle).to_eulerangles_body121<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body131", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body131 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body131<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body131<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body131<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body131<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body131<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body131<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body131(0.5, 0.0, 0.5).to_eulerangles_body131<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body131(1.0, 0.0, 0.0).to_eulerangles_body131<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body131(1.0, PI, 0.0).to_eulerangles_body131<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body131(some_angle, some_angle, some_angle).to_eulerangles_body131<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body212", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body212 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body212<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body212<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body212<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body212<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body212<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body212<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body212(0.5, 0.0, 0.5).to_eulerangles_body212<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body212(1.0, 0.0, 0.0).to_eulerangles_body212<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body212(1.0, PI, 0.0).to_eulerangles_body212<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body212(some_angle, some_angle, some_angle).to_eulerangles_body212<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body232", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body232 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body232<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body232<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body232<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body232<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(-PI/2.0, 1));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body232<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body232<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body232(0.5, 0.0, 0.5).to_eulerangles_body232<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body232(1.0, 0.0, 0.0).to_eulerangles_body232<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body232(1.0, PI, 0.0).to_eulerangles_body232<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body232(some_angle, some_angle, some_angle).to_eulerangles_body232<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body313", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body313 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body313<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body313<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body313<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body313<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body313<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(-PI/2.0, 1));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body313<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body313(0.5, 0.0, 0.5).to_eulerangles_body313<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body313(1.0, 0.0, 0.0).to_eulerangles_body313<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body313(1.0, PI, 0.0).to_eulerangles_body313<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body313(some_angle, some_angle, some_angle).to_eulerangles_body313<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_body323", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body323 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body323<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body323<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body323<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body323<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body323<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body323<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body323(0.5, 0.0, 0.5).to_eulerangles_body323<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body323(1.0, 0.0, 0.0).to_eulerangles_body323<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body323(1.0, PI, 0.0).to_eulerangles_body323<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body323(some_angle, some_angle, some_angle).to_eulerangles_body323<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_space123", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space123 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space123<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space123<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space123<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space123<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space123<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space123<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space123(0.5, PI/2.0, 0.5).to_eulerangles_space123<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space123(1.0, PI/2.0, 0.0).to_eulerangles_space123<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space123(1.0, -PI/2.0, 0.0).to_eulerangles_space123<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space123(some_angle, some_angle, some_angle).to_eulerangles_space123<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}


TEST_CASE("to_eulerangles_space231", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space231 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space231<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space231<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space231<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space231<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space231<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space231<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space231(0.5, PI/2.0, 0.5).to_eulerangles_space231<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space231(1.0, PI/2.0, 0.0).to_eulerangles_space231<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space231(1.0, -PI/2.0, 0.0).to_eulerangles_space231<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space231(some_angle, some_angle, some_angle).to_eulerangles_space231<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}


TEST_CASE("to_eulerangles_space312", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space312 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space312<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space312<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space312<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space312<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space312<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space312<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space312(0.5, PI/2.0, 0.5).to_eulerangles_space312<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space312(1.0, PI/2.0, 0.0).to_eulerangles_space312<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space312(1.0, -PI/2.0, 0.0).to_eulerangles_space312<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space312(some_angle, some_angle, some_angle).to_eulerangles_space312<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_space132", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space132 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space132<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space132<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space132<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space132<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space132<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space132<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space132(0.5, PI/2.0, 0.5).to_eulerangles_space132<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space132(1.0, PI/2.0, 0.0).to_eulerangles_space132<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space132(1.0, -PI/2.0, 0.0).to_eulerangles_space132<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space132(some_angle, some_angle, some_angle).to_eulerangles_space132<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_space213", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space213 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space213<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space213<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space213<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space213<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space213<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space213<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space213(0.5, PI/2.0, 0.5).to_eulerangles_space213<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space213(1.0, PI/2.0, 0.0).to_eulerangles_space213<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space213(1.0, -PI/2.0, 0.0).to_eulerangles_space213<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space213(some_angle, some_angle, some_angle).to_eulerangles_space213<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_space321", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space321 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space321<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space321<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space321<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space321<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space321<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space321<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space321(0.5, PI/2.0, 0.5).to_eulerangles_space321<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space321(1.0, PI/2.0, 0.0).to_eulerangles_space321<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space321(1.0, -PI/2.0, 0.0).to_eulerangles_space321<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space321(some_angle, some_angle, some_angle).to_eulerangles_space321<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_space121", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space121 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space121<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space121<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space121<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space121<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space121<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space121<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space121(0.5, 0.0, 0.5).to_eulerangles_space121<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space121(1.0, 0.0, 0.0).to_eulerangles_space121<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space121(1.0, PI, 0.0).to_eulerangles_space121<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space121(some_angle, some_angle, some_angle).to_eulerangles_space121<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_space131", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space131 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space131<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space131<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space131<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space131<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space131<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(-PI/2.0, 1));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space131<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space131(0.5, 0.0, 0.5).to_eulerangles_space131<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space131(1.0, 0.0, 0.0).to_eulerangles_space131<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space131(1.0, PI, 0.0).to_eulerangles_space131<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space131(some_angle, some_angle, some_angle).to_eulerangles_space131<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_space212", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space212 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space212<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space212<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space212<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space212<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space212<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space212<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(-PI/2.0, 1));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space212(0.5, 0.0, 0.5).to_eulerangles_space212<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space212(1.0, 0.0, 0.0).to_eulerangles_space212<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space212(1.0, PI, 0.0).to_eulerangles_space212<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space212(some_angle, some_angle, some_angle).to_eulerangles_space212<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_space232", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space232 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space232<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space232<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space232<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space232<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space232<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space232<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space232(0.5, 0.0, 0.5).to_eulerangles_space232<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space232(1.0, 0.0, 0.0).to_eulerangles_space232<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space232(1.0, PI, 0.0).to_eulerangles_space232<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space232(some_angle, some_angle, some_angle).to_eulerangles_space232<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_space313", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space313 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space313<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space313<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space313<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space313<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space313<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(-PI/2.0, 1));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(PI/2.0, 1));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space313<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space313(0.5, 0.0, 0.5).to_eulerangles_space313<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space313(1.0, 0.0, 0.0).to_eulerangles_space313<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space313(1.0, PI, 0.0).to_eulerangles_space313<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space313(some_angle, some_angle, some_angle).to_eulerangles_space313<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

TEST_CASE("to_eulerangles_space323", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space323 function.
    */

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space323<double>());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space323<double>());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space323<double>();
    CHECK_THAT(ex1theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex1theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space323<double>();
    CHECK_THAT(ex2theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex2theta3, Catch::Matchers::WithinULP(-PI/2.0, 1));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space323<double>();
    CHECK_THAT(ex3theta1, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex3theta2, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex3theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space323<double>();
    CHECK_THAT(ex4theta1, Catch::Matchers::WithinULP(PI/2.0, 1));
    CHECK_THAT(ex4theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex4theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space323(0.5, 0.0, 0.5).to_eulerangles_space323<double>();
    CHECK_THAT(ex5theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space323(1.0, 0.0, 0.0).to_eulerangles_space323<double>();
    CHECK_THAT(ex5g1theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g1theta2, Catch::Matchers::WithinULP(0.0, 0));
    CHECK_THAT(ex5g1theta3, Catch::Matchers::WithinULP(0.0, 0));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space323(1.0, PI, 0.0).to_eulerangles_space323<double>();
    CHECK_THAT(ex5g2theta1, Catch::Matchers::WithinULP(1.0, 0));
    CHECK_THAT(ex5g2theta2, Catch::Matchers::WithinULP(PI, 0));
    CHECK_THAT(ex5g2theta3, Catch::Matchers::WithinULP(0.0, 0));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space323(some_angle, some_angle, some_angle).to_eulerangles_space323<double>();
    CHECK_THAT(ex6theta1, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta2, Catch::Matchers::WithinULP(some_angle, 1));
    CHECK_THAT(ex6theta3, Catch::Matchers::WithinULP(some_angle, 1));
}

