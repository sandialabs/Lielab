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

#include "pybind_optimize.hpp"
#include "optimize/pybind_ExtremizationMethods.hpp"
#include "optimize/pybind_RootSolvingMethods.hpp"

namespace py = pybind11;

void bind_optimize(py::module& m_optimize)
{
    bind_ExtremizationMethods(m_optimize);
    bind_RootSolvingMethods(m_optimize);
}
