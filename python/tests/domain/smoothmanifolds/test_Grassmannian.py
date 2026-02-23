import pytest

def test_Grassmannian_to_string():
    from lielab.domain import Grassmannian

    xblank = Grassmannian()
    assert (xblank.to_string() == "Grassmannian(0, 0, R)")

    x0 = Grassmannian(0, 0)
    assert (x0.to_string() == "Grassmannian(0, 0, R)")
    x1 = Grassmannian(1, 1)
    assert (x1.to_string() == "Grassmannian(1, 1, R)")
    x10 = Grassmannian(5, 10)
    assert (x10.to_string() == "Grassmannian(5, 10, R)")

def test_Grassmannian_main_initializer():
    from lielab.domain import Grassmannian

    xblank = Grassmannian()
    assert (xblank.get_dimension() == 0)

    x0 = Grassmannian(0, 0)
    assert (x0.get_dimension() == 0)
    x1 = Grassmannian(1, 1)
    assert (x1.get_dimension() == 1)
    x10 = Grassmannian(5, 10)
    assert (x10.get_dimension() == 5)

    with pytest.raises(RuntimeError):
        Grassmannian(6, 5)

def test_Grassmannian_get_dimension():
    from lielab.domain import Grassmannian

    zero = Grassmannian(0, 0)
    one = Grassmannian(1, 1)
    two = Grassmannian(2, 2)
    three = Grassmannian(3, 3)
    four = Grassmannian(4, 4)
    five = Grassmannian(5, 5)
    six = Grassmannian(6, 6)
    seven = Grassmannian(7, 7)
    eight = Grassmannian(8, 8)

    assert (zero.get_dimension() == 0)
    assert (one.get_dimension() == 1)
    assert (two.get_dimension() == 2)
    assert (three.get_dimension() == 3)
    assert (four.get_dimension() == 4)
    assert (five.get_dimension() == 5)
    assert (six.get_dimension() == 6)
    assert (seven.get_dimension() == 7)
    assert (eight.get_dimension() == 8)

def test_Grassmannian_serialize_unserialize():
    """
    Tests the serialize/unserialize operation.
    """

    from lielab.domain import Grassmannian

    x0 = Grassmannian(0, 0)
    x0.unserialize([])
    x0bar = x0.serialize()

    assert (x0bar.size == 0)

    x1 = Grassmannian(1, 2)
    x1.unserialize([1.0, 2.0, 3.0])
    x1bar = x1.serialize()

    assert (x1bar.size == 2)
    assert (x1bar[0] == 1.0)
    assert (x1bar[1] == 2.0)

    x1.unserialize([4.0, 5.0])
    x1bar = x1.serialize()

    assert (x1bar.size == 2)
    assert (x1bar[0] == 4.0)
    assert (x1bar[1] == 5.0)

    x1.unserialize([6.0])
    x1bar = x1.serialize()

    assert (x1bar.size == 2)
    assert (x1bar[0] == 6.0)
    assert (x1bar[1] == 5.0)

    x2 = Grassmannian(2, 4)
    x2.unserialize([1.0, 2.0, 3.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 4)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 3.0)
    assert (x2bar[3] == 0.0)

    x2.unserialize([4.0, 5.0, 6.0, 7.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 4)
    assert (x2bar[0] == 4.0)
    assert (x2bar[1] == 5.0)
    assert (x2bar[2] == 6.0)
    assert (x2bar[3] == 7.0)

    x2.unserialize([8.0, 9.0, 10.0, 11.0, 12.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 4)
    assert (x2bar[0] == 8.0)
    assert (x2bar[1] == 9.0)
    assert (x2bar[2] == 10.0)
    assert (x2bar[3] == 11.0)

def test_Grassmannian_project_point():
    from lielab.domain import Grassmannian
    import numpy as np
    
    pi = np.pi
    axes = np.array([[np.cos(1.0/6.0*pi)], [np.sin(1.0/6.0*pi)]])

    Gr12 = Grassmannian([0.0, 0.0], axes)

    vals1 = np.array([np.cos(1.0/6.0*pi), np.sin(1.0/6.0*pi)])
    proj1 = Gr12.project_point(vals1)

    assert (proj1.size == 2)
    assert (np.abs(proj1[0] - np.cos(1.0/6.0*pi)) < 1e-15)
    assert (np.abs(proj1[1] - np.sin(1.0/6.0*pi)) < 1e-15)

    vals2 = np.array([-np.sin(1.0/6.0*pi), np.cos(1.0/6.0*pi)])
    proj2 = Gr12.project_point(vals2)

    assert (proj2.size == 2)
    assert (np.abs(proj2[0] - 0.0) < 1e-15)
    assert (np.abs(proj2[1] - 0.0) < 1e-15)

    vals3 = np.array([np.cos(1.0/6.0*pi), np.sin(1.0/6.0*pi) - 4.0])
    proj3 = Gr12.project_point(vals3)

    assert (proj3.size == 2)
    assert (np.abs(proj3[0] + np.cos(1.0/6.0*pi)) < 1e-15)
    assert (np.abs(proj3[1] + np.sin(1.0/6.0*pi)) < 1e-15)

def test_Grassmannian_project_vector_onto_tangent_space():
    from lielab.domain import Grassmannian
    import numpy as np
    
    pi = np.pi
    vector = np.array([np.sqrt(2.0)/2.0, np.sqrt(2.0)/2.0])

    axes1 = np.array([[1.0], [0.0]])
    Gr121 = Grassmannian.project([0.0, 0.0], axes1)
    proj1 = Gr121.project_vector_onto_tangent_space(vector)

    assert (proj1.size == 2)
    assert (np.abs(proj1[0] - np.sqrt(2.0)/2.0) < 1e-15)
    assert (np.abs(proj1[1] - 0.0) < 1e-15)

    axes2 = np.array([[0.0], [1.0]])
    Gr122 = Grassmannian.project([0.0, 0.0], axes2)
    proj2 = Gr122.project_vector_onto_tangent_space(vector)

    assert (proj2.size == 2)
    assert (np.abs(proj2[0] - 0.0) < 1e-15)
    assert (np.abs(proj2[1] - np.sqrt(2.0)/2.0) < 1e-15)

    axes3 = np.array([[1.0], [1.0]])
    Gr123 = Grassmannian.project([0.0, 0.0], axes3)
    proj3 = Gr123.project_vector_onto_tangent_space(vector)

    assert (proj3.size == 2)
    assert (np.abs(proj3[0] - np.sqrt(2.0)/2.0) < 1e-15)
    assert (np.abs(proj3[1] - np.sqrt(2.0)/2.0) < 1e-15)

def test_Grassmannian_project_vector_onto_normal_space():
    from lielab.domain import Grassmannian
    import numpy as np
    
    pi = np.pi
    vector = np.array([np.sqrt(2.0)/2.0, np.sqrt(2.0)/2.0])

    axes1 = np.array([[1.0], [0.0]])
    Gr121 = Grassmannian.project([0.0, 0.0], axes1)
    proj1 = Gr121.project_vector_onto_normal_space(vector)

    assert (proj1.size == 2)
    assert (np.abs(proj1[0] - 0.0) < 1e-15)
    assert (np.abs(proj1[1] - np.sqrt(2.0)/2.0) < 1e-15)

    axes2 = np.array([[0.0], [1.0]])
    Gr122 = Grassmannian.project([0.0, 0.0], axes2)
    proj2 = Gr122.project_vector_onto_normal_space(vector)

    assert (proj2.size == 2)
    assert (np.abs(proj2[0] - np.sqrt(2.0)/2.0) < 1e-15)
    assert (np.abs(proj2[1] - 0.0) < 1e-15)

    axes3 = np.array([[1.0], [1.0]])
    Gr123 = Grassmannian.project([0.0, 0.0], axes3)
    proj3 = Gr123.project_vector_onto_normal_space(vector)

    assert (proj3.size == 2)
    assert (np.abs(proj3[0] - 0.0) < 1e-15)
    assert (np.abs(proj3[1] - 0.0) < 1e-15)

def test_Grassmannian_axes_intersection():
    from lielab.domain import Grassmannian
    import numpy as np

    axes1 = np.array([[1.0], [0.0]])
    Gr121 = Grassmannian.project([0.0, 0.0], axes1)
    intersection1 = Gr121.axes_intersection([0.0, 1.0], [0.0, -1.0])

    assert (intersection1.size == 2)
    assert (np.abs(intersection1[0] - 0.0) < 1e-15)
    assert (np.abs(intersection1[1] - 0.0) < 1e-15)

    intersection2 = Gr121.axes_intersection([0.0, 1.0], [1.0, 0.0])

    assert (intersection2.size == 2)
    assert (np.isnan(intersection2[0]))
    assert (np.isnan(intersection2[1]))

    intersection3 = Gr121.axes_intersection([0.0, 1.0], [1.0, -1.0])

    assert (intersection3.size == 2)
    assert (np.abs(intersection3[0] - 1.0) < 1e-15)
    assert (np.abs(intersection3[1] - 0.0) < 1e-15)

    axes2 = np.array([[1.0], [1.0]])
    Gr122 = Grassmannian.project([0.0, 0.0], axes2)
    intersection4 = Gr122.axes_intersection([0.0, 1.0], [0.0, -1.0])

    assert (intersection4.size == 2)
    assert (np.abs(intersection4[0] - 0.0) < 1e-15)
    assert (np.abs(intersection4[1] - 0.0) < 1e-15)

    intersection5 = Gr122.axes_intersection([0.0, 1.0], [1.0, 0.0])

    assert (intersection5.size == 2)
    assert (np.abs(intersection5[0] - 1.0) < 1e-15)
    assert (np.abs(intersection5[1] - 1.0) < 1e-15)

    intersection6 = Gr122.axes_intersection([0.0, 1.0], [1.0, -1.0])

    assert (intersection6.size == 2)
    assert (np.abs(intersection6[0] - 0.5) < 1e-15)
    assert (np.abs(intersection6[1] - 0.5) < 1e-15)

