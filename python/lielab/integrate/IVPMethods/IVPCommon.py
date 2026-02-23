from lielab.cppLielab.integrate import IVPSolution as _IVPSolution
from lielab.cppLielab.integrate import IVPMethod, IVPOptions, IVPStatus
from lielab.cppLielab.integrate import EuclideanIVPSystem as _EuclideanIVPSystem
from lielab.cppLielab.integrate import HomogeneousIVPSystem as _HomogeneousIVPSystem
from lielab.cppLielab.integrate import EuclideanClassicHamiltonianIVPSystem as _EuclideanClassicHamiltonianIVPSystem


class IVPSolution(_IVPSolution):
    def __repr__(self):
        # order = ['messages', 'status', 't', 'y', 'ybar', 'p']
        strout = ''
        strout += '{:>14}: {}\n'.format('message', self.message)
        strout += '{:>14}: {}\n'.format('status', self.status)
        strout += '{:>14}: {}\n'.format('t', 'NumPy Array ' + str(self.t.shape))
        strout += '{:>14}: {}\n'.format('y', 'CompositeManifold (' + str(len(self.y)) + ',)')
        strout += '{:>14}: {}\n'.format('ybar', 'NumPy Array ' + str(self.ybar.shape))
        strout += '{:>14}: {}\n'.format('theta', 'CompositeAlgebra (' + str(len(self.theta)) + ',)')
        strout += '{:>14}: {}\n'.format('thetabar', 'NumPy Array ' + str(self.thetabar.shape))
        return strout

class EuclideanIVPSystem(_EuclideanIVPSystem):
    def __init__(self, vectorfield, event=None):
        super(EuclideanIVPSystem, self).__init__(vectorfield)

        if event is not None:
            self.event = event


class HomogeneousIVPSystem(_HomogeneousIVPSystem):
    def __init__(self, generator,
                 action=None,
                 connection=None,
                 coordinates=None,
                 event=None):
        super(HomogeneousIVPSystem, self).__init__(generator)

        if action is not None:
            self.action = action
        
        if connection is not None:
            self.connection = connection
        
        if coordinates is not None:
            self.coordinates = coordinates
        
        if event is not None:
            self.event = event

class EuclideanClassicHamiltonianIVPSystem(_EuclideanClassicHamiltonianIVPSystem):
    def __init__(self, dVdq, dTdp, event=None):
        super(EuclideanClassicHamiltonianIVPSystem, self).__init__(dVdq, dTdp)

        if event is not None:
            self.event = event