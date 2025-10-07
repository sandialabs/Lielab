from .IVPMethods import *

import numpy as np

def solve_ivp(dynamics, tspan, y0, options = IVPOptions()):
    from lielab.integrate import RungeKuttaFlow
    from lielab.integrate import EuclideanIVPSystem, HomogeneousIVPSystem

    # TODO: Check that len(tspan) > 1 here
    # TODO: Check that tspan is ascending

    if isinstance(dynamics, EuclideanIVPSystem):
        F = RungeKuttaFlow()

        # Preprocess
        if not F.method.can_variable_step:
            F.variable_time_step = False
        
        if options.reltol is not None:
            F.reltol = float(options.reltol)
            F.method.reltol = float(options.reltol)
        
        if options.abstol is not None:
            F.abstol = float(options.abstol)
            F.method.abstol = float(options.abstol)
        
        if options.dt_initial is not None:
            F.dt = float(options.dt_initial)

        if options.dt_min is not None:
            F.dt_min = float(options.dt_min)
        
        if options.dt_max is not None:
            F.dt_max = float(options.dt_max)
        
        if options.small is not None:
            F.small = options.small
        
        if options.large is not None:
            F.large = options.large
        
        if options.pessimist is not None:
            F.pessimist = options.pessimist

        # Run the RK method
        out = F(dynamics, tspan, y0, options)

        # Postprocess
        if (out.status == 0):
            out.message = "Method converged."
            out.success = True
        elif (out.status == -1):
            out.message = "NaNs in vectorfield."
            out.success = False
        elif (out.status == -2):
            out.message = "Infs in vectorfield."
            out.success = False

        return out
    
    elif isinstance(dynamics, HomogeneousIVPSystem):
        xi0 = dynamics.vectorfield(tspan[0], y0)
        rksols = []
        y0seg = [y0]
        for ii in range(len(tspan) - 1):
            def vectorfield_wrapped(t, thetabar, xi0=xi0, y0=y0seg[ii], dynamics=dynamics):
                theta = xi0*0.0 # Awkwardly force a copy
                theta.set_vector(thetabar)
                Theta = dynamics.coordinates(theta)
                y = dynamics.action(Theta, y0)
                xi = dynamics.vectorfield(t, y)
                dy = dynamics.connection(theta, xi)
                return dy.get_vector()

            MuntheKaasZannaDynamics = EuclideanIVPSystem(vectorfield_wrapped)

            event_test = dynamics.event(tspan[0], y0seg[ii])
            if not np.isnan(event_test):
                def event_wrapped(t, thetabar, xi0=xi0, y0=y0seg[ii], dynamics=dynamics):
                    theta = xi0*0.0 # Awkwardly force a copy
                    theta.set_vector(thetabar)
                    Theta = dynamics.coordinates(theta)
                    y = dynamics.action(Theta, y0)
                    return dynamics.event(t, y)

                MuntheKaasZannaDynamics.event = event_wrapped

            # TODO: Feed new dt
            segment = solve_ivp(MuntheKaasZannaDynamics, [tspan[ii], tspan[ii+1]], np.zeros((xi0.get_dimension())), options)

            n_t = segment.t.size
            segment.theta = []
            segment.y = []
            _ybar = np.zeros((n_t, y0.serialize().size))

            for jj in range(len(segment.t)):
                thetaj = 0.0*xi0
                thetaj.set_vector(segment.ybar[jj, :])
                segment.theta += [thetaj]
                Thetaj = dynamics.coordinates(thetaj)
                yj = dynamics.action(Thetaj, y0seg[ii])
                segment.y += [yj]
                _ybar[jj, :] = yj.serialize()

            segment.ybar = _ybar
            y0seg += [segment.y[-1]]
            rksols += [segment]

        out = rksols[0] # TODO: Make a copy?

        # Concatenate each segment into a single curve
        for ii in range(len(tspan) - 2):
            # TODO: The overlapping point is removed. Should it be kept in?
            out.t = np.hstack((out.t, rksols[ii+1].t[1:]))
            out.theta += rksols[ii+1].theta[1:]
            out.thetabar = np.vstack((out.thetabar, rksols[ii+1].thetabar[1:, :]))
            out.y += rksols[ii+1].y[1:]
            out.ybar = np.vstack((out.ybar, rksols[ii+1].ybar[1:, :]))

        return out

    else:
        raise NotImplementedError()

    
