from lielab.testing import *

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
