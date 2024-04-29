from lielab import ALGO_STATUS
from lielab.domain import halie, hmlie
from lielab.topos import IntegralCurve
from lielab.cppLielab.topos import Flow as _Flow
from lielab.cppLielab.topos import MuntheKaas as _MuntheKaas
from lielab.cppLielab.topos import TimeStepper as _TimeStepper
from lielab.cppLielab.topos import TSOutput
import numpy as np

# Import hnewton as a python object
from lielab.optim import hnewton

class MuntheKaas(_MuntheKaas):
    """
    MuntheKaas integrator that maximizes time spent in C++
    for performance.
    """
    
    def __call__(self, vectorfield, y0, t0, dt):
        self.init(t0, y0, dt, vectorfield(t0, y0))

        while (self.algo_status == ALGO_STATUS.OK):
            self.step_0()
            self.set_dy(vectorfield(self.next_t, self.next_y))
            self.step_1()

        return self.postprocess()


class CustomMuntheKaas(_TimeStepper):
    """
    MuntheKaas integrator with customizable actions.
    """

    def __init__(self, *args, **kwargs):
        """
        Instantiates a new MuntheKaas object.
        """

        super(CustomMuntheKaas, self).__init__(*args, **kwargs)

        from lielab.functions import left_product
        from lielab.topos import exp, dexpinv
        # self.RKT = RKTYPE::RKTYPE_EXPLICIT; # TODO:

        self.action = left_product
        self.phi = exp
        self.dphiinv = dexpinv

        self._t0 = None
        self._y0 = None
        self._dt = None

        self.next_t = None
        self.next_y = None

        self._dy = None
        self._U = None
        self._KK = None

    def init(self, t0, y0, dt, dy0):
        """
        Initializes the MuntheKaas process.
        """
        
        super(CustomMuntheKaas, self).init()

        self._t0 = t0
        self._y0 = y0
        self._dt = dt
        self._KK = [halie()]*self.n
        self._KK[0] = dy0
        self._dy = dy0
        self._U = 0*dy0

        if (self.iterations >= self.n - 1):
            self.algo_status = ALGO_STATUS.FINISHED
        
    def step_0(self):
        """
        Advances the solution to current iteration and returns pair
        (next_t, next_y) to be evaluated by the vectorfield.
        """

        self._U *= 0

        for jj in range(self.iterations + 1):
            self._U += self._dt*self.A[self.iterations+1, jj]*self._KK[jj]
        
        self.next_y = self.action(self.phi(self._U), self._y0)
        self.next_t = self._t0 + self._dt*self.C[self.iterations + 1]
    
    def set_dy(self, dy):
        """
        Sets the current vectorfield value. Do this after step_0()
        and before step_1().
        """

        self._dy = dy

    def step_1(self):
        """
        Evaluates the $ d phi^{-1} u ( xi(s,y)) $ map.
        """

        self._KK[self.iterations + 1] = self.dphiinv(self._U, self._dy, self.n - 1)

        super(CustomMuntheKaas, self).step()

        if (self.iterations == self.n - 1):
            self.algo_status = ALGO_STATUS.FINISHED
    
    def postprocess(self):
        Ulow = 0*self._KK[0]
        Uhigh = 0*self._KK[0]

        for ii in range(self.n):
            Ulow += self._dt*self.B[ii]*self._KK[ii]

        low = self.action(self.phi(Ulow), self._y0)

        if self.variable_step:
            for ii in range(self.n):
                Uhigh += self._dt*self.Bhat[ii]*self._KK[ii]
            
            high = self.action(self.phi(Uhigh), self._y0)
            error = np.linalg.norm(Uhigh.get_vector() - Ulow.get_vector())
            return TSOutput(low, high, error)
        
        return TSOutput(low)
    
    def __call__(self, vectorfield, y0, t0, dt):
        self.init(t0, y0, dt, vectorfield(t0, y0))

        while (self.algo_status == ALGO_STATUS.OK):
            self.step_0()
            self.set_dy(vectorfield(self.next_t, self.next_y))
            self.step_1()
        
        return self.postprocess()


class Flow(_Flow):
    search = hnewton()
    stepper = MuntheKaas()

    def __call__(self, vectorfield, tspan, y0, *args):
        from lielab.domain import rn, RN, halie, hmlie

        # Check if the solution we're running has an event
        has_event = False
        if len(args) > 0:
            has_event = True
            event = args[0]
        
        # Check the type of inputs
        dy0 = vectorfield(tspan[0], y0)

        if isinstance(dy0, np.ndarray) or isinstance(dy0, list):
            # NumPy and list type inputs, we need to wrap these in
            # sophus objects
            _y0 = hmlie([RN(y0)])
            def _vectorfield(_t, M):
                _RN = M.space[0]
                _states = _RN.serialize()
                gradient = vectorfield(_t, _states)
                dx = rn(gradient)
                return halie([dx])
            
            if has_event:
                def _event(_t, M):
                    _RN = M.space[0]
                    _states = _RN.serialize()
                    return event(_t, _states)
        else:
            # Assume these are sophus type inputs and
            # don't need to be wrapped
            _y0 = y0
            _vectorfield = vectorfield
            if has_event:
                _event = event
        
        # Initialize the Flow object
        self.init(tspan, _y0)

        # Main loop. Run until algorithm says it's done
        while (self.algo_status == ALGO_STATUS.OK):
            accepted = False

            # This is the RK loop. It's the error control loop that will
            # change step size until we have an acceptable step
            while (not accepted):
                next = self.stepper(_vectorfield, self._ynext, self._out.t[self.iterations], self._dt)
                if has_event:
                    if _event(self._out.t[self.iterations], self._ynext) > 0 and _event(self._out.t[self.iterations] + self._dt, next.high) <= 0:
                        self.search.lower = halie([rn([self.tol])])
                        self.search.upper = halie([rn([self._dt])])

                        def sfun(sdt):
                            _sdt = sdt.space[0](0)
                            return _event(self._out.t[self.iterations] + _sdt, self.stepper(_vectorfield, self._ynext, self._out.t[self.iterations], _sdt).high)

                        temp_dt = self.search(sfun, halie([rn([self._dt/2])]))
                        self._dt = temp_dt.space[0](0)
                        next = self.stepper(_vectorfield, self._ynext, self._out.t[self.iterations], self._dt)

                accepted = self.step0(next)
            
            # We have an acceptable step at this point. Step the solution forward by one
            self.step()

            # If there's an event, check to see if conditions are met
            # This will stop integration if an event returns <= 0
            if has_event:
                self.stepE(_event(self._out.t[self.iterations], self._ynext))
        
        self.postprocess()
        
        return IntegralCurve(self._out)
