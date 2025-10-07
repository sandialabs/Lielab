#include <Lielab.hpp>
#include <iostream>
#include <numbers>

#include <catch2/catch_all.hpp>

#include "../../test_utils.hpp"

TEST_CASE("from_SO3", "[domain]")
{
    /*!
    * Tests the from_SO3 function
    * TODO: update to new SU serialization
    */
    
    constexpr double PI = std::numbers::pi_v<double>;

    // Build 90 degree rotations in x, y, and z
    Lielab::domain::SO rx = Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(0,3));
    Lielab::domain::SO ry = Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(1,3));
    Lielab::domain::SO rz = Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(2,3));

    // Test 90 degree x rotation
    Lielab::domain::SU _qx = Lielab::domain::SU::from_SO3(rx);
    std::array<double, 4> qx = _qx.to_quaternion();

    CHECK(std::abs(qx[0] - std::sqrt(2.0)/2.0) <= TOL_FINE);
    CHECK(std::abs(qx[1] - std::sqrt(2.0)/2.0) <= TOL_FINE);
    CHECK(std::abs(qx[2] - 0.0) <= TOL_FINE);
    CHECK(std::abs(qx[3] - 0.0) <= TOL_FINE);

    // Test 90 degree y rotation
    Lielab::domain::SU _qy = Lielab::domain::SU::from_SO3(ry);
    std::array<double, 4> qy = _qy.to_quaternion();

    CHECK(std::abs(qy[0] - std::sqrt(2.0)/2.0) <= TOL_FINE);
    CHECK(std::abs(qy[1] - 0.0) <= TOL_FINE);
    CHECK(std::abs(qy[2] - std::sqrt(2.0)/2.0) <= TOL_FINE);
    CHECK(std::abs(qy[3] - 0.0) <= TOL_FINE);

    // Test 90 degree x rotation
    Lielab::domain::SU _qz = Lielab::domain::SU::from_SO3(rz);
    std::array<double, 4> qz = _qz.to_quaternion();

    CHECK(std::abs(qz[0] - std::sqrt(2.0)/2.0) <= TOL_FINE);
    CHECK(std::abs(qz[1] - 0.0) <= TOL_FINE);
    CHECK(std::abs(qz[2] - 0.0) <= TOL_FINE);
    CHECK(std::abs(qz[3] - std::sqrt(2.0)/2.0) <= TOL_FINE);

}

TEST_CASE("from_quaternion", "[domain]")
{
    /*!
    * Tests quaternions against well-known identities.
    */

    // Lielab::domain::SU q1 = Lielab::domain::SU::from_quaternion();
    Lielab::domain::SU qm1 = Lielab::domain::SU::from_quaternion(-1.0, 0.0, 0.0, 0.0);
    Lielab::domain::SU qi = Lielab::domain::SU::from_quaternion(0.0, 1.0, 0.0, 0.0);
    Lielab::domain::SU qj = Lielab::domain::SU::from_quaternion(0.0, 0.0, 1.0, 0.0);
    Lielab::domain::SU qk = Lielab::domain::SU::from_quaternion(0.0, 0.0, 0.0, 1.0);

    // std::vector<Lielab::domain::SU> elements;
    // elements.push_back(qi);
    // elements.push_back(qj);
    // elements.push_back(qk);

    // is_group<Lielab::domain::SU>(elements, q1);

    // Hamilton's identities
    // i^2 = j^2 = k^2 = -1
    assert_domain(qi*qi, qm1);
    assert_domain(qj*qj, qm1);
    assert_domain(qk*qk, qm1);

    // ij = -ji = -k
    assert_domain(qi*qj, (qj*qi).inverse());
    assert_domain(qi*qj, qk.inverse());

    // jk = -kj = -i
    assert_domain(qj*qk, (qk*qj).inverse());
    assert_domain(qj*qk, qi.inverse());

    // ki = -ik = -j
    assert_domain(qk*qi, (qi*qk).inverse());
    assert_domain(qk*qi, qj.inverse());
}
