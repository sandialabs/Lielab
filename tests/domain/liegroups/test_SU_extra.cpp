#include <Lielab/domain/liegroups/SU.hpp>
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

    using Lielab::domain::SU;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;

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
    CHECK(check_almost_equal_nulp(CompositeGroup{qi*qi}, CompositeGroup{qm1}, 1, true));
    CHECK(check_almost_equal_nulp(CompositeGroup{qj*qj}, CompositeGroup{qm1}, 1, true));
    CHECK(check_almost_equal_nulp(CompositeGroup{qk*qk}, CompositeGroup{qm1}, 1, true));

    // ij = -ji = -k
    CHECK(check_almost_equal_nulp(CompositeGroup{qi*qj}, CompositeGroup{(qj*qi).inverse()}, 1, true));
    CHECK(check_almost_equal_nulp(CompositeGroup{qi*qj}, CompositeGroup{qk.inverse()}, 1, true));

    // jk = -kj = -i
    CHECK(check_almost_equal_nulp(CompositeGroup{qj*qk}, CompositeGroup{(qk*qj).inverse()}, 1, true));
    CHECK(check_almost_equal_nulp(CompositeGroup{qj*qk}, CompositeGroup{qi.inverse()}, 1, true));

    // ki = -ik = -j
    CHECK(check_almost_equal_nulp(CompositeGroup{qk*qi}, CompositeGroup{(qi*qk).inverse()}, 1, true));
    CHECK(check_almost_equal_nulp(CompositeGroup{qk*qi}, CompositeGroup{qj.inverse()}, 1, true));
}
