from lielab.cppLielab.integrate import IVPOptions
from lielab.cppLielab.integrate import EuclideanIVPSystem as _EuclideanIVPSystem
from lielab.cppLielab.integrate import HomogeneousIVPSystem as _HomogeneousIVPSystem

class EuclideanIVPSystem(_EuclideanIVPSystem):
    def __init__(self, vectorfield, event=None):
        super(EuclideanIVPSystem, self).__init__(vectorfield)

        if event is not None:
            self.event = event


class HomogeneousIVPSystem(_HomogeneousIVPSystem):
    def __init__(self, vectorfield,
                 action=None,
                 connection=None,
                 coordinates=None,
                 event=None):
        super(HomogeneousIVPSystem, self).__init__(vectorfield)

        if action is not None:
            self.action = action
        
        if connection is not None:
            self.connection = connection
        
        if coordinates is not None:
            self.coordinates = coordinates
        
        if event is not None:
            self.event = event
