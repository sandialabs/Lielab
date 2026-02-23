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

#include "pybind_RootSolvingMethods.hpp"
#include "RootSolving/pybind_RootSolvingCommon.hpp"
#include "RootSolving/pybind_NewtonRootSearch.hpp"
#include "RootSolving/pybind_HybridRootSearch.hpp"

namespace py = pybind11;

void bind_RootSolvingMethods(py::module& m_optimize)
{
    bind_RootSolvingCommon(m_optimize);
    bind_NewtonRootSearch(m_optimize);
    bind_HybridRootSearch(m_optimize);
}
