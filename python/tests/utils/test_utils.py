from lielab.testing import *

def optimfun(x):
    return x[0]*x[0]*x[0] - 2*x[0] - 5

def test_opt_golden():
    from lielab.utils import opt_golden

    search = opt_golden()

    search.lower = [0.0]
    search.upper = [2.0]

    search.init()

    assert abs(search.lower[0] - 0.0) < TOL_FINE
    assert abs(search.upper[0] - 2.0) < TOL_FINE
    assert abs(search._A[0] - 0.0) < TOL_FINE
    assert abs(search._B[0] - 2.0) < TOL_FINE
    assert abs(search._X1[0] - 0.763932022500210) < TOL_FINE
    assert abs(search._X2[0] - 1.236067977499790) < TOL_FINE

    search._f1 = optimfun(search._X1)
    search.num_objective_evals += 1
    search._f2 = optimfun(search._X2)
    search.num_objective_evals += 1
    search.step()

    assert abs(search.lower[0] - 0.0) < TOL_FINE
    assert abs(search.upper[0] - 2.0) < TOL_FINE
    assert abs(search._A[0] - 0.0) < TOL_FINE
    assert abs(search._B[0] - 1.236067977499790) < TOL_FINE
    assert abs(search._X1[0] - 0.472135954999579) < TOL_FINE
    assert abs(search._X2[0] - 0.763932022500210) < TOL_FINE

    search._f1 = optimfun(search._X1)
    search.num_objective_evals += 1
    search._f2 = optimfun(search._X2)
    search.num_objective_evals += 1
    search.step()

    assert abs(search.lower[0] - 0.0) < TOL_FINE
    assert abs(search.upper[0] - 2.0) < TOL_FINE
    assert abs(search._A[0] - 0.472135954999579) < TOL_FINE
    assert abs(search._B[0] - 1.236067977499790) < TOL_FINE
    assert abs(search._X1[0] - 0.763932022500210) < TOL_FINE
    assert abs(search._X2[0] - 0.944271909999159) < TOL_FINE

    search.init()

    assert search.iterations == 0
    assert search.num_objective_evals == 0
    assert search.num_jacobian_evals == 0
    assert search.num_hessian_evals == 0

    assert abs(search.lower[0] - 0.0) < TOL_FINE
    assert abs(search.upper[0] - 2.0) < TOL_FINE
    assert abs(search._A[0] - 0.0) < TOL_FINE
    assert abs(search._B[0] - 2.0) < TOL_FINE
    assert abs(search._X1[0] - 0.763932022500210) < TOL_FINE
    assert abs(search._X2[0] - 1.236067977499790) < TOL_FINE

    xopt = search(optimfun)

    assert abs(xopt[0] - 0.816496623231839) < TOL_FINE
    assert search.iterations == 31
    # assert search.success == True
    # assert search.algo_status == ALGO_STATUS.FINISHED

    search.init()
    search.max_iterations = 5

    xopt = search(optimfun)

    # assert search.success == False
    assert search.iterations == 5
    # assert search.algo_status == ALGO_STATUS.MAXITER


def optfun2(x):
    return x**2 - 1

def test_newton():
    from lielab.utils import newton

    y01 = 0.5

    search = newton()

    yopt = search(optfun2, y01)

    assert abs(yopt - 1.0) <= TOL_FINE

# def test_search_linearx():
#     from lielab.optim import search_linearx
#     search = search_linearx()
#     fun = lambda x: x - 4.5
#     x = 5.0

#     search.lower = 4.0
#     search.upper = 6.0

#     search.init(x)

#     assert abs(search.lower - 4.0) < TOL_FINE
#     assert abs(search.upper - 6.0) < TOL_FINE

#     x = search.step(x, fun(x))

#     assert abs(x - 5.000006) < TOL_FINE
#     assert search.k == 3

#     x = search.step(x, fun(x))

#     assert abs(x - 4.5) < TOL_FINE
#     assert search.k == 4

#     x = search.step(x, fun(x))

#     assert abs(x - 4.5) < TOL_FINE
#     assert search.k == 8

#     x = 5.0

#     x = search(fun, x)

#     assert abs(x - 4.5) < TOL_FINE
#     assert search.k == 8
