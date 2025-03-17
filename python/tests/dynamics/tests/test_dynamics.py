import numpy as np

from lielab.testing import *

# def test_vfex1():
#     """
#     Tests vfex1.
#     """

#     from lielab.domain import SO, CompositeManifold
#     from lielab.dynamics import vfex1

#     vf = vfex1()
#     y0 = SO(3)
#     M0 = CompositeManifold([y0])

#     out = vf(2.0, M0)

#     dy = out.space[0]

#     abs(dy(0,0) - 0.0) < TOL_FINE
#     abs(dy(0,1) - 2.0) < TOL_FINE
#     abs(dy(0,2) - 1.0) < TOL_FINE
#     abs(dy(1,0) + 2.0) < TOL_FINE
#     abs(dy(1,1) - 0.0) < TOL_FINE
#     abs(dy(1,2) + 4.0) < TOL_FINE
#     abs(dy(2,0) + 1.0) < TOL_FINE
#     abs(dy(2,1) - 4.0) < TOL_FINE
#     abs(dy(2,2) - 0.0) < TOL_FINE
