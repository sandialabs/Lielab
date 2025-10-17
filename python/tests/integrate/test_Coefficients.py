
def test_Coefficients():
    from lielab.integrate import Coefficients, get_butcher_tableau
    import math
    import numpy as np

    for method in Coefficients:
        [A, b, bhat, c, e, order, stages, variable, implicit] = get_butcher_tableau(method)

        assert (A.shape[0] == stages)
        assert (b.shape[0] == stages)

        if (method != Coefficients.RKV65e):
            err = np.abs(np.sum(b) - 1.0)
            assert (err <= math.ulp(1))

        if (variable):
            assert (bhat.shape[0] == stages)

            if (method != Coefficients.RKV87e):
                err = np.abs(np.sum(bhat) - 1.0)
                assert (err <= math.ulp(1))
        
        assert (c.shape[0] == stages)
        assert (e.shape[0] == stages)
