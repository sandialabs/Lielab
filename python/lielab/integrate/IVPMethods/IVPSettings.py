from lielab.cppLielab.integrate import IVPOptions

from lielab.functions import left_Lie_group_action, dexpinv, exp

import numpy as np

# class IVPOptions(_IVPOptions):
#     def __init__(self, *args, **kwargs):
#         super(IVPOptions, self).__init__(*args, **kwargs)

class EuclideanIVPSystem(object):
    def __init__(self, vectorfield, event=None):

        if event is not None:
            self.__dict__.update({'event': event})
        else:
            self.__dict__.update({'event': lambda _t, _y: np.nan})
        self.__dict__.update({'vectorfield': vectorfield})

class HomogeneousIVPSystem(object):
    def __init__(self, vectorfield,
                 action=left_Lie_group_action,
                 connection=dexpinv,
                 coordinates=exp,
                 event=None):
        
        self.__dict__.update({'action': action})
        self.__dict__.update({'connection': connection})
        self.__dict__.update({'coordinates': coordinates})
        if event is not None:
            self.__dict__.update({'event': event})
        else:
            self.__dict__.update({'event': lambda _t, _y: np.nan})
        self.__dict__.update({'vectorfield': vectorfield})
