from lielab.cppLielab.optimize import ExtremizationMethod, ExtremizationOptions
from lielab.cppLielab.optimize import EuclideanExtremizationSystem as _EuclideanExtremizationSystem

class EuclideanExtremizationSystem(_EuclideanExtremizationSystem):
    def __init__(self, objective, jacobian=None, hessian=None):
        super(EuclideanExtremizationSystem, self).__init__(objective)

        if jacobian is not None:
            self.jacobian = jacobian
        
        if hessian is not None:
            self.hessian = hessian
