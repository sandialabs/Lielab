from lielab.cppLielab.integrate import RungeKuttaStatus

from .Coefficients import Coefficients
from lielab.cppLielab.integrate import RungeKutta as _RungeKutta

from lielab.cppLielab.integrate import RungeKuttaFlowStatus
from lielab.cppLielab.integrate import RungeKuttaFlow as RungeKuttaFlowCPP

from lielab.integrate import ODESolution
from lielab.utils import newton


import numpy as np

class RungeKutta(_RungeKutta):
    """
    """
    def __init__(self, scheme=Coefficients.RKV87r):
        """
        Instantiates a new RungeKutta object.
        """
        super(RungeKutta, self).__init__(scheme)

    def __call__(self, vectorfield, y0, t0, dt):
        self.status = 0
        self.message = ''

        dy = vectorfield(t0, y0)
        status = super(RungeKutta, self).init(t0, dt, dy)

        while (status.value > 0):
            if (status == RungeKuttaStatus.DO_STEP0):
                status = super(RungeKutta, self).step_0()
            elif (status == RungeKuttaStatus.DO_STEP1):
                next_y = y0 + self.next_theta
                dy = vectorfield(self.next_t, next_y)
                if np.any(np.isnan(dy)):
                    self.status = -1
                    self.message = 'NaNs in vectorfield.'
                    return y0
                if np.any(np.isinf(dy)):
                    self.status = -2
                    self.message = 'Infs in vectorfield.'
                    return y0

                status = super(RungeKutta, self).step_1(dy)

        status = super(RungeKutta, self).postprocess()
        next_y = y0 + self.next_theta
        if (status == RungeKuttaStatus.ESTIMATE_ERROR):
            status = self.estimate_error(y0 + self.next_theta, y0 + self.next_theta2)
        
        return next_y

class RungeKuttaFlow(RungeKuttaFlowCPP):
    search = newton()
    method = RungeKutta()

    def __call__(self, dynamics, tspan, y0, options):
        import numpy as np
        from lielab.domain import RN, CompositeManifold

        # Check if the solution we're running has an event
        event_val = dynamics.event(tspan[0], y0)
        self.has_event = False
        if not np.isnan(event_val):
            self.has_event = True

        # Wrap data
        _y0 = y0

        if isinstance(_y0, list):
            _y0 = np.array(_y0)

        # Check the type of inputs
        dy0 = dynamics.vectorfield(tspan[0], _y0)

        # Initialize the Flow object
        status = self.init(tspan, _y0)

        # Main loop. Run until algorithm says it's done
        while (status.value > 0):
            if (status == RungeKuttaFlowStatus.DO_STEP0):
                ynext = self.method(dynamics.vectorfield, self._ycurrent, self._tcurrent, self.dt)
                if self.method.status != 0:
                    break

                if self.has_event:
                    self.event_current = dynamics.event(self._tcurrent, self._ycurrent)
                    self.event_next = dynamics.event(self._tcurrent + self.dt, ynext)
                    if (self.event_current >= 0) and (self.event_next <= 0):
                        # Event crossed
                        self.search.lower = np.finfo(np.float64).eps
                        self.search.upper = self.dt
                        def fun(_dt):
                            _ynext = self.method(dynamics.vectorfield, self._ycurrent, self._tcurrent, _dt)
                            return dynamics.event(self._tcurrent + _dt, _ynext)
                        xopt = self.search(fun, self.dt/2.0)
                        self.dt = xopt
                        ynext = self.method(dynamics.vectorfield, self._ycurrent, self._tcurrent, self.dt)
                        self.event_next = dynamics.event(self._tcurrent + self.dt, ynext)

                status = self.step0(ynext, self.method.error_estimate)
            elif (status == RungeKuttaFlowStatus.DO_STEP1):
                status = self.step1()

        self.postprocess()

        out = ODESolution(self._out)
        out.status = self.method.status

        return out
