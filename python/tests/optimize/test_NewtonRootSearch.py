
def test_NewtonRootSearch_requires_jacobian():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.5]

    def objective(x):
        return [np.pow(x[0], 2.0) - 1.0]
    
    system = EuclideanRootSystem(objective)

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == False)

def test_NewtonRootSearch_nans_objective():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.5]

    def objective(x):
        return [np.pow(x[0], 2.0)*np.nan - 1.0]
    
    def jacobian(x):
        return [[2.0*x[0]]]
    
    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == False)


def test_NewtonRootSearch_infs_objective():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.5]

    def objective(x):
        return [np.pow(x[0], 2.0)*np.inf - 1.0]
    
    def jacobian(x):
        return [[2.0*x[0]]]
    
    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == False)

def test_NewtonRootSearch_nans_jacobian():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.5]

    def objective(x):
        return [np.pow(x[0], 2.0) - 1.0]
    
    def jacobian(x):
        return [[2.0*x[0]*np.nan]]
    
    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == False)


def test_NewtonRootSearch_infs_jacobian():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.5]

    def objective(x):
        return [np.pow(x[0], 2.0) - 1.0]
    
    def jacobian(x):
        return [[2.0*x[0]*np.inf]]
    
    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == False)


def test_NewtonRootSearch_1D_free():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.5]

    def objective(x):
        return [np.pow(x[0], 2.0) - 1.0]

    def jacobian(x):
        return [[2.0*x[0]]]
    
    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == True)
    assert (np.abs(xopt[0] - 1.0) <= options.tol)


def test_NewtonRootSearch_1D_lower():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.5]

    def objective(x):
        return [np.pow(x[0], 2.0) - 1.0]

    def jacobian(x):
        return [[2.0*x[0]]]
    
    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian
    system.lower_bound = [1.1]

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == False)
    assert (np.abs(xopt[0] - 1.1) <= options.tol)


def test_NewtonRootSearch_1D_upper():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.5]

    def objective(x):
        return [np.pow(x[0], 2.0) - 1.0]
    
    def jacobian(x):
        return [[2.0*x[0]]]
    
    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian
    system.upper_bound = [0.9]

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == False)
    assert (np.abs(xopt[0] - 0.9) <= options.tol)


def test_NewtonRootSearch_1D_bounded():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.5]

    def objective(x):
        return [np.pow(x[0], 2.0) - 1.0]
    
    def jacobian(x):
        return [[2.0*x[0]]]
    
    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian
    system.lower_bound = [0.9]
    system.upper_bound = [1.1]

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == True)
    assert (np.abs(xopt[0] - 1.0) <= options.tol)


def test_NewtonRootSearch_2D_free():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.0, 0.0]

    A = np.zeros((2, 2))
    A[0, 0] = 1.0
    A[0, 1] = 2.0
    A[1, 0] = -3.0
    A[1, 1] = 4.0

    B = np.array([0.5, 0.5])

    def objective(x):
        return np.dot(A, x) + B
    
    def jacobian(x):
        return A
    
    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == True)
    assert (np.abs(xopt[0] + 0.1) <= options.tol)
    assert (np.abs(xopt[1] + 0.2) <= options.tol)
    assert (solver.iteration <= 1)


def test_NewtonRootSearch_2D_lower():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.0, 0.0]

    A = np.zeros((2, 2))
    A[0, 0] = 1.0
    A[0, 1] = 2.0
    A[1, 0] = -3.0
    A[1, 1] = 4.0

    B = np.array([0.5, 0.5])

    def objective(x):
        return np.dot(A, x) + B
    
    def jacobian(x):
        return A

    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian
    system.lower_bound = [0.0]

    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == False)
    assert (np.abs(xopt[0] + 0.0) <= options.tol)
    assert (np.abs(xopt[1] + 0.2) <= options.tol)
    assert (solver.iteration <= 5)

def test_NewtonRootSearch_3D_underdetermined():
    from lielab.optimize import NewtonRootSearch, RootSolvingOptions, RootSolvingMethod, EuclideanRootSystem
    import numpy as np

    x_guess = [0.0, 0.0, 0.0]

    A = np.zeros((2, 3))
    A[0, 0] = 1.0
    A[0, 1] = 2.0
    A[0, 2] = 4.0
    A[1, 0] = -3.0
    A[1, 1] = 4.0
    A[1, 2] = -5.0

    B = np.array([0.5, 0.5])

    def objective(x):
        return np.dot(A, x) + B
    
    def jacobian(x):
        return A
    
    system = EuclideanRootSystem(objective)
    system.jacobian = jacobian
    
    options = RootSolvingOptions()

    solver = NewtonRootSearch()
    xopt = solver(system, x_guess, options)

    assert (solver.success == True)
    objfinal = objective(xopt)
    assert (np.abs(objfinal[0]) <= options.tol)
    assert (np.abs(objfinal[1]) <= options.tol)
    assert (solver.iteration <= 1)

