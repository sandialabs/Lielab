
def test_RungeKuttaCoefficients():
    from lielab.integrate import RungeKuttaCoefficients, get_butcher_tableau
    import math
    import numpy as np

    for method in RungeKuttaCoefficients:
        [A, b, bhat, c, e, order, stages, variable, implicit] = get_butcher_tableau(method)

        for ii in range(stages):
            for jj in range(stages):
                assert not np.isnan(A[ii][jj])

            assert not np.isnan(b[ii])

            if variable:
                assert not np.isnan(bhat[ii])
                assert not np.isnan(e[ii])
            
            assert not np.isnan(c[ii])

        if (method != RungeKuttaCoefficients.RKV65e):
            err = np.abs(np.sum(b) - 1.0)
            # assert (err <= math.ulp(1))

        if (variable):
            if (method != RungeKuttaCoefficients.RKV87e):
                err = np.abs(np.sum(bhat) - 1.0)
                # assert (err <= math.ulp(1))
