def test_sign():
    from lielab.utils import sign
    import numpy as np

    inf = np.inf

    assert (sign(0.25) == sign(0.5))
    assert (sign(-0.25) == sign(-0.5))
    assert (sign(0.25) != sign(-0.5))
    assert (sign(-0.25) != sign(0.5))
    
    assert (sign(0.5) != sign(0.0))
    assert (sign(-0.5) != sign(0.0))
    assert (sign(0.0) == sign(0.0))

    assert (sign(inf) == sign(0.5))
    assert (sign(-inf) == sign(-0.5))
    assert (sign(inf) != sign(-0.5))
    assert (sign(-inf) != sign(0.5))
    assert (sign(inf) == sign(inf))
    assert (sign(-inf) == sign(-inf))
