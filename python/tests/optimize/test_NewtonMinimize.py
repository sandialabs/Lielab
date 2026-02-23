import pytest

def test_NewtonMinimize_requires_jacobian():
    from lielab.optimize import NewtonMinimize, ExtremizationOptions, EuclideanExtremizationSystem
    import numpy as np

    x_guess = [2.0, 2.0]

    def objective(x):
        return (x[0] - 1.0)**2 + (x[1] - 1.0)**2
    
    system = EuclideanExtremizationSystem(objective)

    options = ExtremizationOptions()

    solver = NewtonMinimize()

    with pytest.raises(RuntimeError):
        xopt = solver(system, x_guess, options)

def test_NewtonMinimize_requires_hessian():
    from lielab.optimize import NewtonMinimize, ExtremizationOptions, EuclideanExtremizationSystem
    import numpy as np

    x_guess = [2.0, 2.0]

    def objective(x):
        return x[0]**2 + (x[1] - 0.5)**3 + x[0]*x[1]
    
    def jacobian(x):
        return np.array([2.0*(x[0] - 1.0), 2.0*(x[1] - 1.0)])
    
    system = EuclideanExtremizationSystem(objective)
    system.jacobian = jacobian

    options = ExtremizationOptions()

    solver = NewtonMinimize()
    with pytest.raises(RuntimeError):
        xopt = solver(system, x_guess, options)

def test_NewtonMinimize_size_hessian():
    from lielab.optimize import NewtonMinimize, ExtremizationOptions, EuclideanExtremizationSystem
    import numpy as np

    x_guess = [2.0, 2.0]

    def objective(x):
        return x[0]**2 + (x[1] - 0.5)**3 + x[0]*x[1]
    
    def jacobian(x):
        return np.array([2.0*(x[0] - 1.0), 2.0*(x[1] - 1.0)])
    
    def hessian(x):
        out = np.zeros((3,3))
        out[0,0] = 2.0
        out[0,1] = 0.0
        out[1,0] = 0.0
        out[1,1] = 2.0
        return out
    
    system = EuclideanExtremizationSystem(objective)
    system.jacobian = jacobian
    system.hessian = hessian

    options = ExtremizationOptions()

    solver = NewtonMinimize()
    with pytest.raises(RuntimeError):
        xopt = solver(system, x_guess, options)

def test_NewtonMinimize_simple():
    from lielab.optimize import NewtonMinimize, ExtremizationOptions, EuclideanExtremizationSystem
    import numpy as np

    x_guess = [2.0, 2.0]

    def objective(x):
        return x[0]**2 + (x[1] - 0.5)**3 + x[0]*x[1]
    
    def jacobian(x):
        return np.array([2.0*(x[0] - 1.0), 2.0*(x[1] - 1.0)])
    
    def hessian(x):
        out = np.zeros((2,2))
        out[0,0] = 2.0
        out[0,1] = 0.0
        out[1,0] = 0.0
        out[1,1] = 2.0
        return out
    
    system = EuclideanExtremizationSystem(objective)
    system.jacobian = jacobian
    system.hessian = hessian

    options = ExtremizationOptions()

    solver = NewtonMinimize()
    xopt = solver(system, x_guess, options)

    assert solver.success == True
    assert np.abs(xopt[0] - 1.0) < options.abstol
    assert np.abs(xopt[0] - 1.0) < options.abstol
