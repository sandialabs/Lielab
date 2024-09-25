#include <Lielab.hpp>
#include <iostream>

#include "test_utils.hpp"

TEST_CASE("from_quaternion", "[domain]")
{
    /*!
    * Tests quaternions against well-known identities.
    */

    Lielab::domain::SU q1 = Lielab::domain::SU::from_quaternion();
    Lielab::domain::SU qm1 = Lielab::domain::SU::from_quaternion(-1.0, 0.0, 0.0, 0.0);
    Lielab::domain::SU qi = Lielab::domain::SU::from_quaternion(0.0, 1.0, 0.0, 0.0);
    Lielab::domain::SU qj = Lielab::domain::SU::from_quaternion(0.0, 0.0, 1.0, 0.0);
    Lielab::domain::SU qk = Lielab::domain::SU::from_quaternion(0.0, 0.0, 0.0, 1.0);

    std::vector<Lielab::domain::SU> elements;
    elements.push_back(qi);
    elements.push_back(qj);
    elements.push_back(qk);

    is_group<Lielab::domain::SU>(elements, q1);

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
