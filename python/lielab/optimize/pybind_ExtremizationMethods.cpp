#include <Lielab.hpp>
#include <string>
#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/native_enum.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "pybind_ExtremizationMethods.hpp"
#include "Extremization/pybind_ExtremizationCommon.hpp"
#include "Extremization/pybind_LineSearch.hpp"
#include "Extremization/pybind_GoldenMinimize.hpp"
#include "Extremization/pybind_NewtonMinimize.hpp"

namespace py = pybind11;

void bind_ExtremizationMethods(py::module& m_optimize)
{
    bind_ExtremizationCommon(m_optimize);
    bind_LineSearch(m_optimize);
    bind_GoldenMinimize(m_optimize);
    bind_NewtonMinimize(m_optimize);
}
