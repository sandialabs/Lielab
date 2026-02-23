
def test_GoldenMinimize_fail_bounds():
    from lielab.optimize import GoldenMinimize, EuclideanExtremizationSystem, ExtremizationOptions, ExtremizationMethod
    import numpy as np

    def optimfun(x):
        return x[0]*x[0]*x[0] - 2*x[0] - 5

    system = EuclideanExtremizationSystem(optimfun)

    options = ExtremizationOptions()

    solver = GoldenMinimize()
    xopt = solver(system, options)

    assert (solver.success == False)

def test_GoldenMinimize_1D():
    from lielab.optimize import GoldenMinimize, EuclideanExtremizationSystem, ExtremizationOptions, ExtremizationMethod
    import numpy as np

    def optimfun(x):
        return x[0]*x[0]*x[0] - 2*x[0] - 5

    system = EuclideanExtremizationSystem(optimfun)
    system.lower_bound = np.array([0.0])
    system.upper_bound = np.array([2.0])

    options = ExtremizationOptions()

    solver = GoldenMinimize()
    xopt = solver(system, options)

    assert (solver.success == True)
    assert (np.abs(xopt[0] - np.sqrt(2.0/3.0)) < options.reltol)

def test_GoldenMinimize_1D_right_nans():
    from lielab.optimize import GoldenMinimize, EuclideanExtremizationSystem, ExtremizationOptions, ExtremizationMethod
    import numpy as np

    def optimfun(x):
        if x[0] > 0.8:
            return np.nan
        return 1 - x[0]

    system = EuclideanExtremizationSystem(optimfun)
    system.lower_bound = np.array([0.0])
    system.upper_bound = np.array([2.0])

    options = ExtremizationOptions()

    solver = GoldenMinimize()
    xopt = solver(system, options)

    assert (solver.success == True)
    assert (np.abs(xopt[0] - 0.8) < options.reltol)

def test_GoldenMinimize_1D_left_nans():
    from lielab.optimize import GoldenMinimize, EuclideanExtremizationSystem, ExtremizationOptions, ExtremizationMethod
    import numpy as np

    def optimfun(x):
        if x[0] < 1.2:
            return np.nan
        return x[0]

    system = EuclideanExtremizationSystem(optimfun)
    system.lower_bound = np.array([0.0])
    system.upper_bound = np.array([2.0])

    options = ExtremizationOptions()

    solver = GoldenMinimize()
    xopt = solver(system, options)

    assert (solver.success == True)
    assert (np.abs(xopt[0] - 1.2) < options.reltol)
