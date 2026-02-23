#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

#include <complex>

TEST_CASE("left_Lie_group_action", "[functions]")
{
    /*!
    * Tests the left_Lie_group_action function.
    */
    
    using namespace Lielab::domain;
    using namespace Lielab::functions;

    const std::complex<double> j(0.0, 1.0);

    // Wrong sizes
    CHECK_THROWS(left_Lie_group_action(CompositeGroup({SO(3), SU(3)}), CompositeManifold({SO(3)})));
    CHECK_THROWS(left_Lie_group_action(CompositeGroup({SO(3)}), CompositeManifold({SO(3), SU(3)})));

    // Unknown action
    CHECK_THROWS(left_Lie_group_action(CompositeGroup({SO(3)}), CompositeManifold({SU(3)})));

    // CN x CN
    const CN g_CN = CN::from_complex_vector({1.0 + 2.0*j, 3.0 + 4.0*j});
    const CN y_CN = CN(2);

    const CompositeManifold test_CN_CN = left_Lie_group_action(CompositeGroup({g_CN}), CompositeManifold({y_CN}));
    REQUIRE(test_CN_CN.point.size() == 1);
    const CN test_CN_out = test_CN_CN[0];
    const Eigen::VectorXcd test_CN_outbar = test_CN_out.to_complex_vector();
    REQUIRE(test_CN_outbar.size() == 2);
    CHECK(test_CN_outbar(0) == std::complex<double>(1.0, 2.0));
    CHECK(test_CN_outbar(1) == std::complex<double>(3.0, 4.0));

    // GLR x GLR
    const CompositeManifold test_GLR_GLR = left_Lie_group_action(CompositeGroup({GLR(2)}), CompositeManifold({GLR(2)}));
    REQUIRE(test_GLR_GLR.point.size() == 1);
    const GLR test_GLR_out = test_GLR_GLR[0];
    CHECK(test_GLR_out.get_shape() == 2);

    // GLC x GLC
    const CompositeManifold test_GLC_GLC = left_Lie_group_action(CompositeGroup({GLC(2)}), CompositeManifold({GLC(2)}));
    REQUIRE(test_GLC_GLC.point.size() == 1);
    const GLC test_GLC_out = test_GLC_GLC[0];
    CHECK(test_GLC_out.get_shape() == 2);

    // RN x RN
    const RN g_RN = RN::from_vector({1.0, 2.0, 3.0, 4.0});
    const RN y_RN = RN(4);

    const CompositeManifold test_RN_RN = left_Lie_group_action(CompositeGroup({g_RN}), CompositeManifold({y_RN}));
    REQUIRE(test_RN_RN.point.size() == 1);
    const RN test_RN_out = test_RN_RN[0];
    const Eigen::VectorXd test_RN_outbar = test_RN_out.serialize();
    REQUIRE(test_RN_outbar.size() == 4);
    CHECK(test_RN_outbar(0) == 1.0);
    CHECK(test_RN_outbar(1) == 2.0);
    CHECK(test_RN_outbar(2) == 3.0);
    CHECK(test_RN_outbar(3) == 4.0);

    // SE x SE
    const CompositeManifold test_SE_SE = left_Lie_group_action(CompositeGroup({SE(2)}), CompositeManifold({SE(2)}));
    REQUIRE(test_SE_SE.point.size() == 1);
    const SE test_SE_out = test_SE_SE[0];
    CHECK(test_SE_out.get_shape() == 3);

    // SO x SO
    const CompositeManifold test_SO_SO = left_Lie_group_action(CompositeGroup({SO(2)}), CompositeManifold({SO(2)}));
    REQUIRE(test_SO_SO.point.size() == 1);
    const SO test_SO_out = test_SO_SO[0];
    CHECK(test_SO_out.get_shape() == 2);

    // SP x SP
    const CompositeManifold test_SP_SP = left_Lie_group_action(CompositeGroup({SP(2)}), CompositeManifold({SP(2)}));
    REQUIRE(test_SP_SP.point.size() == 1);
    const SP test_SP_out = test_SP_SP[0];
    CHECK(test_SP_out.get_shape() == 2);

    // SU x SU
    const CompositeManifold test_SU_SU = left_Lie_group_action(CompositeGroup({SU(2)}), CompositeManifold({SU(2)}));
    REQUIRE(test_SU_SU.point.size() == 1);
    const SU test_SU_out = test_SU_SU[0];
    CHECK(test_SU_out.get_shape() == 2);
}

TEST_CASE("right_Lie_group_action", "[functions]")
{
    /*!
    * Tests the right_Lie_group_action function.
    */
    
    using namespace Lielab::domain;
    using namespace Lielab::functions;

    const std::complex<double> j(0.0, 1.0);

    // Wrong sizes
    CHECK_THROWS(right_Lie_group_action(CompositeGroup({SO(3), SU(3)}), CompositeManifold({SO(3)})));
    CHECK_THROWS(right_Lie_group_action(CompositeGroup({SO(3)}), CompositeManifold({SO(3), SU(3)})));

    // Unknown action
    CHECK_THROWS(right_Lie_group_action(CompositeGroup({SO(3)}), CompositeManifold({SU(3)})));

    // CN x CN
    const CN g_CN = CN::from_complex_vector({1.0 + 2.0*j, 3.0 + 4.0*j});
    const CN y_CN = CN(2);

    const CompositeManifold test_CN_CN = right_Lie_group_action(CompositeGroup({g_CN}), CompositeManifold({y_CN}));
    REQUIRE(test_CN_CN.point.size() == 1);
    const CN test_CN_out = test_CN_CN[0];
    const Eigen::VectorXcd test_CN_outbar = test_CN_out.to_complex_vector();
    REQUIRE(test_CN_outbar.size() == 2);
    CHECK(test_CN_outbar(0) == std::complex<double>(1.0, 2.0));
    CHECK(test_CN_outbar(1) == std::complex<double>(3.0, 4.0));

    // GLR x GLR
    const CompositeManifold test_GLR_GLR = right_Lie_group_action(CompositeGroup({GLR(2)}), CompositeManifold({GLR(2)}));
    REQUIRE(test_GLR_GLR.point.size() == 1);
    const GLR test_GLR_out = test_GLR_GLR[0];
    CHECK(test_GLR_out.get_shape() == 2);

    // GLC x GLC
    const CompositeManifold test_GLC_GLC = right_Lie_group_action(CompositeGroup({GLC(2)}), CompositeManifold({GLC(2)}));
    REQUIRE(test_GLC_GLC.point.size() == 1);
    const GLC test_GLC_out = test_GLC_GLC[0];
    CHECK(test_GLC_out.get_shape() == 2);

    // RN x RN
    const RN g_RN = RN::from_vector({1.0, 2.0, 3.0, 4.0});
    const RN y_RN = RN(4);

    const CompositeManifold test_RN_RN = right_Lie_group_action(CompositeGroup({g_RN}), CompositeManifold({y_RN}));
    REQUIRE(test_RN_RN.point.size() == 1);
    const RN test_RN_out = test_RN_RN[0];
    const Eigen::VectorXd test_RN_outbar = test_RN_out.serialize();
    REQUIRE(test_RN_outbar.size() == 4);
    CHECK(test_RN_outbar(0) == 1.0);
    CHECK(test_RN_outbar(1) == 2.0);
    CHECK(test_RN_outbar(2) == 3.0);
    CHECK(test_RN_outbar(3) == 4.0);

    // SE x SE
    const CompositeManifold test_SE_SE = right_Lie_group_action(CompositeGroup({SE(2)}), CompositeManifold({SE(2)}));
    REQUIRE(test_SE_SE.point.size() == 1);
    const SE test_SE_out = test_SE_SE[0];
    CHECK(test_SE_out.get_shape() == 3);

    // SO x SO
    const CompositeManifold test_SO_SO = right_Lie_group_action(CompositeGroup({SO(2)}), CompositeManifold({SO(2)}));
    REQUIRE(test_SO_SO.point.size() == 1);
    const SO test_SO_out = test_SO_SO[0];
    CHECK(test_SO_out.get_shape() == 2);

    // SP x SP
    const CompositeManifold test_SP_SP = right_Lie_group_action(CompositeGroup({SP(2)}), CompositeManifold({SP(2)}));
    REQUIRE(test_SP_SP.point.size() == 1);
    const SP test_SP_out = test_SP_SP[0];
    CHECK(test_SP_out.get_shape() == 2);

    // SU x SU
    const CompositeManifold test_SU_SU = right_Lie_group_action(CompositeGroup({SU(2)}), CompositeManifold({SU(2)}));
    REQUIRE(test_SU_SU.point.size() == 1);
    const SU test_SU_out = test_SU_SU[0];
    CHECK(test_SU_out.get_shape() == 2);
}
