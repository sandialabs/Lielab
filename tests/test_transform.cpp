#include <cmath>
#include "lielab.hpp"

constexpr double PI = lielab::constants::PI<double>;
const lielab::domain::SO DCMId(3);
const lielab::domain::SO DCMrotx(lielab::functions::exp(PI/2.0*lielab::domain::so::basis(0,3)));
const lielab::domain::SO DCMroty(lielab::functions::exp(PI/2.0*lielab::domain::so::basis(1,3)));
const lielab::domain::SO DCMrotz(lielab::functions::exp(PI/2.0*lielab::domain::so::basis(2,3)));
const double some_angle = PI/2.0*5.0/7.0;

// TEST_CASE("dcm_to_quaternion", "[transform]")
// {
//     /*!
//     * Tests the dcm_to_quaternion function
//     * TODO: update to new SU serialization
//     */

//     // Build 90 degree rotations in x, y, and z
//     lielab::domain::SO rx = lielab::functions::exp(PI/2.0*lielab::domain::so::basis(0,3));
//     lielab::domain::SO ry = lielab::functions::exp(PI/2.0*lielab::domain::so::basis(1,3));
//     lielab::domain::SO rz = lielab::functions::exp(PI/2.0*lielab::domain::so::basis(2,3));

//     // Test 90 degree x rotation
//     lielab::domain::SU _qx = lielab::transform::dcm_to_quaternion(rx);
//     Eigen::VectorXd qx = _qx.serialize();

//     CHECK(std::abs(qx(0) - std::sqrt(2.0)/2.0) <= TOL_FINE);
//     CHECK(std::abs(qx(1) - std::sqrt(2.0)/2.0) <= TOL_FINE);
//     CHECK(std::abs(qx(2) - 0.0) <= TOL_FINE);
//     CHECK(std::abs(qx(3) - 0.0) <= TOL_FINE);

//     // Test 90 degree y rotation
//     lielab::domain::SU _qy = lielab::transform::dcm_to_quaternion(ry);
//     Eigen::VectorXd qy = _qy.serialize();

//     CHECK(std::abs(qy(0) - std::sqrt(2.0)/2.0) <= TOL_FINE);
//     CHECK(std::abs(qy(1) - 0.0) <= TOL_FINE);
//     CHECK(std::abs(qy(2) - std::sqrt(2.0)/2.0) <= TOL_FINE);
//     CHECK(std::abs(qy(3) - 0.0) <= TOL_FINE);

//     // Test 90 degree x rotation
//     lielab::domain::SU _qz = lielab::transform::dcm_to_quaternion(rz);
//     Eigen::VectorXd qz = _qz.serialize();

//     CHECK(std::abs(qz(0) - std::sqrt(2.0)/2.0) <= TOL_FINE);
//     CHECK(std::abs(qz(1) - 0.0) <= TOL_FINE);
//     CHECK(std::abs(qz(2) - 0.0) <= TOL_FINE);
//     CHECK(std::abs(qz(3) - std::sqrt(2.0)/2.0) <= TOL_FINE);

// }

TEST_CASE("dcm_to_eanglebody123", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody123 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody123<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody123<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody123<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody123<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - PI/2.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody123<double>(lielab::transform::eanglebody123_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody123<double>(lielab::transform::eanglebody123_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody231", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody231 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody231<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody231<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody231<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody231<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody231<double>(lielab::transform::eanglebody231_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody231<double>(lielab::transform::eanglebody231_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody312", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody312 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody312<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody312<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody312<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody312<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody312<double>(lielab::transform::eanglebody312_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody312<double>(lielab::transform::eanglebody312_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody132", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody132 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody132<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody132<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody132<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody132<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody132<double>(lielab::transform::eanglebody132_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody132<double>(lielab::transform::eanglebody132_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody213", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody213 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody213<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody213<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody213<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody213<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - PI/2.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody213<double>(lielab::transform::eanglebody213_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody213<double>(lielab::transform::eanglebody213_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody321", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody321 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody321<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody321<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody321<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody321<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody321<double>(lielab::transform::eanglebody321_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody321<double>(lielab::transform::eanglebody321_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody121", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody121 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody121<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody121<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody121<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody121<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 + PI/2.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody121<double>(lielab::transform::eanglebody121_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody121<double>(lielab::transform::eanglebody121_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody131", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody131 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody131<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody131<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody131<double>(DCMroty);
    CHECK(std::abs(ex3theta1 + PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody131<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody131<double>(lielab::transform::eanglebody131_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody131<double>(lielab::transform::eanglebody131_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody212", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody212 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody212<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody212<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody212<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody212<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 + PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - PI/2.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody212<double>(lielab::transform::eanglebody212_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody212<double>(lielab::transform::eanglebody212_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody232", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody232 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody232<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody232<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 + PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody232<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody232<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody232<double>(lielab::transform::eanglebody232_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody232<double>(lielab::transform::eanglebody232_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody313", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody313 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody313<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody313<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody313<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 + PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody313<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody313<double>(lielab::transform::eanglebody313_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody313<double>(lielab::transform::eanglebody313_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglebody323", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglebody323 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglebody323<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglebody323<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 + PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglebody323<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglebody323<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglebody323<double>(lielab::transform::eanglebody323_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglebody323<double>(lielab::transform::eanglebody323_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglespace123", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace123 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace123<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace123<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace123<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace123<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - PI/2.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace123<double>(lielab::transform::eanglespace123_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace123<double>(lielab::transform::eanglespace123_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}


TEST_CASE("dcm_to_eanglespace231", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace231 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace231<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace231<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace231<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace231<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace231<double>(lielab::transform::eanglespace231_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace231<double>(lielab::transform::eanglespace231_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}


TEST_CASE("dcm_to_eanglespace312", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace312 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace312<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace312<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace312<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace312<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace312<double>(lielab::transform::eanglespace312_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace312<double>(lielab::transform::eanglespace312_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglespace132", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglesspace132 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace132<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace132<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace132<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace132<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace132<double>(lielab::transform::eanglespace132_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace132<double>(lielab::transform::eanglespace132_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglespace213", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace213 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace213<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace213<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace213<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace213<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - PI/2.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace213<double>(lielab::transform::eanglespace213_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace213<double>(lielab::transform::eanglespace213_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglespace321", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace321 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace321<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace321<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace321<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace321<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace321<double>(lielab::transform::eanglespace321_to_dcm(0.5, PI/2.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - PI/2.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace321<double>(lielab::transform::eanglespace321_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglespace121", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace121 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace121<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace121<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace121<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace121<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 + PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - PI/2.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace121<double>(lielab::transform::eanglespace121_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace121<double>(lielab::transform::eanglespace121_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglespace131", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace131 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace131<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace131<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace131<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 + PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace131<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace131<double>(lielab::transform::eanglespace131_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace131<double>(lielab::transform::eanglespace131_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglespace212", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace212 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace212<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace212<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace212<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace212<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 + PI/2.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace212<double>(lielab::transform::eanglespace212_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace212<double>(lielab::transform::eanglespace212_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglespace232", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace232 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace232<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace232<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 + PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace232<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace232<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace232<double>(lielab::transform::eanglespace232_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace232<double>(lielab::transform::eanglespace232_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglespace313", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace313 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace313<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace313<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace313<double>(DCMroty);
    CHECK(std::abs(ex3theta1 + PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace313<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace313<double>(lielab::transform::eanglespace313_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace313<double>(lielab::transform::eanglespace313_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}

TEST_CASE("dcm_to_eanglespace323", "[transform]")
{
    /*!
    * Tests the dcm_to_eanglespace323 function.
    */

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = lielab::transform::dcm_to_eanglespace323<double>(DCMId);
    CHECK(ex1theta1 == 0.0);
    CHECK(ex1theta2 == 0.0);
    CHECK(ex1theta3 == 0.0);

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = lielab::transform::dcm_to_eanglespace323<double>(DCMrotx);
    CHECK(std::abs(ex2theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex2theta3 + PI/2.0) <= TOL_FINE);

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = lielab::transform::dcm_to_eanglespace323<double>(DCMroty);
    CHECK(std::abs(ex3theta1 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta2 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex3theta3 - 0.0) <= TOL_FINE);

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = lielab::transform::dcm_to_eanglespace323<double>(DCMrotz);
    CHECK(std::abs(ex4theta1 - PI/2.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta2 - 0.0) <= TOL_FINE);
    CHECK(std::abs(ex4theta3 - 0.0) <= TOL_FINE);

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = lielab::transform::dcm_to_eanglespace323<double>(lielab::transform::eanglespace323_to_dcm(0.5, 0.0, 0.5));
    CHECK(std::abs(ex5theta1 - 1.0) <= TOL_FINE);
    CHECK(std::abs(ex5theta2 - 0.0) <= TOL_FINE);
    CHECK(ex5theta3 == 0.0);

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = lielab::transform::dcm_to_eanglespace323<double>(lielab::transform::eanglespace323_to_dcm(some_angle, some_angle, some_angle));
    CHECK(ex6theta1 == some_angle);
    CHECK(ex6theta2 == some_angle);
    CHECK(ex6theta3 == some_angle);
}


TEST_CASE("eanglebody123_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody123_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody123_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody123_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody123_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody123_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody123_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody231_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody231_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody231_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody231_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody231_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody231_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody231_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody312_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody312_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody312_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody312_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody312_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody312_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody312_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody132_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody132_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody132_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody132_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody132_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody132_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody132_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody213_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody213_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody213_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody213_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody213_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody213_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody213_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody321_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody321_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody321_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody321_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody321_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody321_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody321_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody121_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody121_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody121_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody121_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody121_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody121_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody121_to_dcm<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody121_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody131_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody131_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody131_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody131_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody131_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody131_to_dcm<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody131_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody131_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody212_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody212_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody212_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody212_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody212_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody212_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody212_to_dcm<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody212_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody232_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody232_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody232_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody232_to_dcm<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody232_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody232_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody232_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody232_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody313_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody313_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody313_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody313_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody313_to_dcm<double>(-PI/2.0, -PI/2.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody313_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody313_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody313_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglebody323_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglebody323_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglebody323_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglebody323_to_dcm<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglebody323_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody323_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglebody323_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglebody323_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace123_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace123_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace123_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace123_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace123_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace123_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace123_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace231_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace231_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace231_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace231_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace231_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace231_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace231_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace312_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace231_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace312_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace312_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace312_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace312_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace312_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace132_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace132_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace132_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace132_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace132_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace132_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace132_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace213_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace213_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace213_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace213_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace213_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace213_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace213_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace321_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace321_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace321_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace321_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace321_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace321_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace321_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace121_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace121_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace121_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace121_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace121_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace121_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace121_to_dcm<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace121_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace131_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace131_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace131_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace131_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace131_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace131_to_dcm<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace131_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace131_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace212_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace212_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace212_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace212_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace212_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace212_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace212_to_dcm<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace212_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace232_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace232_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace232_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace232_to_dcm<double>(PI/2.0, -PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace232_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace232_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace232_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace232_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace313_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace313_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace313_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace313_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace313_to_dcm<double>(-PI/2.0, PI/2.0, PI/2.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace313_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace313_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace313_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("eanglespace323_to_dcm", "[transform]")
{
    /*!
    * Tests the eanglespace323_to_dcm function.
    */

    // Identity
    lielab::domain::SO dcm = lielab::transform::eanglespace323_to_dcm<double>(0.0, 0.0, 0.0);
    assert_domain(DCMId, dcm);

    // Rotate 90 degrees by x-axis
    dcm = lielab::transform::eanglespace323_to_dcm<double>(PI/2.0, PI/2.0, -PI/2.0);
    assert_domain(DCMrotx, dcm);

    // Rotate 90 degrees by y-axis
    dcm = lielab::transform::eanglespace323_to_dcm<double>(0.0, PI/2.0, 0.0);
    assert_domain(DCMroty, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace323_to_dcm<double>(PI/2.0, 0.0, 0.0);
    assert_domain(DCMrotz, dcm);

    // Rotate 90 degrees by z-axis
    dcm = lielab::transform::eanglespace323_to_dcm<double>(0.0, 0.0, PI/2.0);
    assert_domain(DCMrotz, dcm);

    // Rotate by a random angle
    dcm = lielab::transform::eanglespace323_to_dcm<double>(some_angle, some_angle, some_angle);
    CHECK(std::abs(dcm.get_ados_representation().determinant() - 1) < TOL_FINE);
    CHECK(std::abs((dcm.get_ados_representation()*dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE);
}

TEST_CASE("quaternion_to_dcm", "[transform]")
{
    /*!
    * Tests the quaternion_to_dcm function
    */

    // Build 90 degree rotations in x, y, and z
    // TODO: simplify with basis()
    lielab::domain::su u(2), v(2), w(2);
    Eigen::Vector3d xx, yy, zz;
    xx << 0.0, 0.0, 0.5;
    yy << 0.0, 0.5, 0.0;
    zz << 0.5, 0.0, 0.0;
    u.set_vector(xx);
    v.set_vector(yy);
    w.set_vector(zz);
    lielab::domain::SU qx = lielab::functions::exp(PI/2.0*u);
    lielab::domain::SU qy = lielab::functions::exp(PI/2.0*v);
    lielab::domain::SU qz = lielab::functions::exp(PI/2.0*w);

    // Test 90 degree x rotation
    lielab::domain::SO rx = lielab::transform::quaternion_to_dcm(qx);

    assert_domain(rx, DCMrotx);

    // Test 90 degree y rotation
    lielab::domain::SO ry = lielab::transform::quaternion_to_dcm(qy);

    assert_domain(ry, DCMroty);

    // Test 90 degree z rotation
    lielab::domain::SO rz = lielab::transform::quaternion_to_dcm(qz);

    assert_domain(rz, DCMrotz);
}
