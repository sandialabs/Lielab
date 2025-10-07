from ..cppLielab.utils import opt_golden as _opt_golden
from ..cppLielab.utils import newton as _newton

import numpy as np

class opt_golden(_opt_golden):
    def __call__(self, function):
        status = self.init()

        while (status != 0 and status != 1):
            self._f1 = function(self._X1)
            self.num_objective_evals += 1

            self._f2 = function(self._X2)
            self.num_objective_evals += 1
            status = self.step()
        
        self._f1 = function(self._X1)
        self.num_objective_evals += 1

        self._f2 = function(self._X2)
        self.num_objective_evals += 1

        if (self._f1 < self._f2):
            self.val_objective = self._f1
            self._X = self._X1
        else:
            self.val_objective = self._f2
            self._X = self._X2

        return self._X

class newton(_newton):
    def __call__(self, fun, m0):
        status = self.init(m0, fun(m0))

        while (status != 0 and status != 1):
            status = self.step0()
            self.fnext = fun(self.mnext)
            status = self.step1()
            self.fnext = fun(self.mnext)

        self.m = self.mnext

        return self.m

# class search_linearx(_search_linearx):
#     def __call__(self, function, x0):
#         self.init(x0)
#         x = x0

#         while (self.algo_status == ALGO_STATUS.OK):
#             x = self.step(x, function(x))

#         return x
