#include <Lielab.hpp>
#include <string>
#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "pybind_special.hpp"

namespace py = pybind11;

void bind_special(py::module &m_utils)
{
    m_utils.def("bernoulli", &Lielab::utils::bernoulli);
    m_utils.def("sinc", &Lielab::utils::sinc);
    m_utils.def("sign", &Lielab::utils::sign);
}
