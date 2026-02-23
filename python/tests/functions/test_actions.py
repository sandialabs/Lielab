def complex(a,b):
    return a + b*1j

def test_left_Lie_group_action():
    from lielab.domain import CN, GLC, GLR, RN, SE, SO, SP, SU, CompositeGroup, CompositeManifold
    from lielab.functions import left_Lie_group_action
    import pytest

    # Wrong sizes
    with pytest.raises(Exception):
        left_Lie_group_action(CompositeGroup([SO(3), SU(3)]), CompositeManifold([SO(3)]))
    
    with pytest.raises(Exception):
        left_Lie_group_action(CompositeGroup([SO(3)]), CompositeManifold([SO(3), SU(3)]))

    # Unknown action
    with pytest.raises(Exception):
        left_Lie_group_action(CompositeGroup([SO(3)]), CompositeManifold([SU(3)]))

    # CN x CN
    g_CN = CN.from_complex_vector([1.0 + 2.0j, 3.0 + 4.0j])
    y_CN = CN(2)

    test_CN_CN = left_Lie_group_action(CompositeGroup([g_CN]), CompositeManifold([y_CN]))
    assert len(test_CN_CN.point) == 1
    test_CN_out = test_CN_CN[0]
    test_CN_outbar = test_CN_out.to_complex_vector()
    assert test_CN_outbar.size == 2
    assert test_CN_outbar[0] == complex(1.0, 2.0)
    assert test_CN_outbar[1] == complex(3.0, 4.0)

    # GLR x GLR
    test_GLR_GLR = left_Lie_group_action(CompositeGroup([GLR(2)]), CompositeManifold([GLR(2)]))
    assert len(test_GLR_GLR.point) == 1
    test_GLR_out = test_GLR_GLR[0]
    assert test_GLR_out.get_shape() == 2

    # GLC x GLC
    test_GLC_GLC = left_Lie_group_action(CompositeGroup([GLC(2)]), CompositeManifold([GLC(2)]))
    assert len(test_GLC_GLC.point) == 1
    test_GLC_out = test_GLC_GLC[0]
    assert test_GLC_out.get_shape() == 2

    # RN x RN
    g_RN = RN.from_vector([1.0, 2.0, 3.0, 4.0])
    y_RN = RN(4)

    test_RN_RN = left_Lie_group_action(CompositeGroup([g_RN]), CompositeManifold([y_RN]))
    assert len(test_RN_RN.point) == 1
    test_RN_out = test_RN_RN[0]
    test_RN_outbar = test_RN_out.serialize()
    assert test_RN_outbar.size == 4
    assert test_RN_outbar[0] == 1.0
    assert test_RN_outbar[1] == 2.0
    assert test_RN_outbar[2] == 3.0
    assert test_RN_outbar[3] == 4.0

    # SE x SE
    test_SE_SE = left_Lie_group_action(CompositeGroup([SE(2)]), CompositeManifold([SE(2)]))
    assert len(test_SE_SE.point) == 1
    test_SE_out = test_SE_SE[0]
    assert test_SE_out.get_shape() == 3

    # SO x SO
    test_SO_SO = left_Lie_group_action(CompositeGroup([SO(2)]), CompositeManifold([SO(2)]))
    assert len(test_SO_SO.point) == 1
    test_SO_out = test_SO_SO[0]
    assert test_SO_out.get_shape() == 2

    # SP x SP
    test_SP_SP = left_Lie_group_action(CompositeGroup([SP(2)]), CompositeManifold([SP(2)]))
    assert len(test_SP_SP.point) == 1
    test_SP_out = test_SP_SP[0]
    assert test_SP_out.get_shape() == 2

    # SU x SU
    test_SU_SU = left_Lie_group_action(CompositeGroup([SU(2)]), CompositeManifold([SU(2)]))
    assert len(test_SU_SU.point) == 1
    test_SU_out = test_SU_SU[0]
    assert test_SU_out.get_shape() == 2


def test_right_Lie_group_action():
    from lielab.domain import CN, GLC, GLR, RN, SE, SO, SP, SU, CompositeGroup, CompositeManifold
    from lielab.functions import right_Lie_group_action
    import pytest

    # Wrong sizes
    with pytest.raises(Exception):
        right_Lie_group_action(CompositeGroup([SO(3), SU(3)]), CompositeManifold([SO(3)]))
    
    with pytest.raises(Exception):
        right_Lie_group_action(CompositeGroup([SO(3)]), CompositeManifold([SO(3), SU(3)]))

    # Unknown action
    with pytest.raises(Exception):
        right_Lie_group_action(CompositeGroup([SO(3)]), CompositeManifold([SU(3)]))

    # CN x CN
    g_CN = CN.from_complex_vector([1.0 + 2.0j, 3.0 + 4.0j])
    y_CN = CN(2)

    test_CN_CN = right_Lie_group_action(CompositeGroup([g_CN]), CompositeManifold([y_CN]))
    assert len(test_CN_CN.point) == 1
    test_CN_out = test_CN_CN[0]
    test_CN_outbar = test_CN_out.to_complex_vector()
    assert test_CN_outbar.size == 2
    assert test_CN_outbar[0] == complex(1.0, 2.0)
    assert test_CN_outbar[1] == complex(3.0, 4.0)

    # GLR x GLR
    test_GLR_GLR = right_Lie_group_action(CompositeGroup([GLR(2)]), CompositeManifold([GLR(2)]))
    assert len(test_GLR_GLR.point) == 1
    test_GLR_out = test_GLR_GLR[0]
    assert test_GLR_out.get_shape() == 2

    # GLC x GLC
    test_GLC_GLC = right_Lie_group_action(CompositeGroup([GLC(2)]), CompositeManifold([GLC(2)]))
    assert len(test_GLC_GLC.point) == 1
    test_GLC_out = test_GLC_GLC[0]
    assert test_GLC_out.get_shape() == 2

    # RN x RN
    g_RN = RN.from_vector([1.0, 2.0, 3.0, 4.0])
    y_RN = RN(4)

    test_RN_RN = right_Lie_group_action(CompositeGroup([g_RN]), CompositeManifold([y_RN]))
    assert len(test_RN_RN.point) == 1
    test_RN_out = test_RN_RN[0]
    test_RN_outbar = test_RN_out.serialize()
    assert test_RN_outbar.size == 4
    assert test_RN_outbar[0] == 1.0
    assert test_RN_outbar[1] == 2.0
    assert test_RN_outbar[2] == 3.0
    assert test_RN_outbar[3] == 4.0

    # SE x SE
    test_SE_SE = right_Lie_group_action(CompositeGroup([SE(2)]), CompositeManifold([SE(2)]))
    assert len(test_SE_SE.point) == 1
    test_SE_out = test_SE_SE[0]
    assert test_SE_out.get_shape() == 3

    # SO x SO
    test_SO_SO = right_Lie_group_action(CompositeGroup([SO(2)]), CompositeManifold([SO(2)]))
    assert len(test_SO_SO.point) == 1
    test_SO_out = test_SO_SO[0]
    assert test_SO_out.get_shape() == 2

    # SP x SP
    test_SP_SP = right_Lie_group_action(CompositeGroup([SP(2)]), CompositeManifold([SP(2)]))
    assert len(test_SP_SP.point) == 1
    test_SP_out = test_SP_SP[0]
    assert test_SP_out.get_shape() == 2

    # SU x SU
    test_SU_SU = right_Lie_group_action(CompositeGroup([SU(2)]), CompositeManifold([SU(2)]))
    assert len(test_SU_SU.point) == 1
    test_SU_out = test_SU_SU[0]
    assert test_SU_out.get_shape() == 2
