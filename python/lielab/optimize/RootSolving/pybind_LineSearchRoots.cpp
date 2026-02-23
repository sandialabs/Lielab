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

#include "pybind_LineSearchRoots.hpp"

namespace py = pybind11;

void bind_LineSearchRoots(py::module& m_optimize)
{
    // auto Lielab_optimize_LineSearchRoots = py::class_<Lielab::optimize::LineSearchRoots>(m_optimize, "LineSearchRoots");
    // Lielab_optimize_LineSearchRoots.def_readwrite("iteration", &Lielab::optimize::LineSearchRoots::iteration);
    // Lielab_optimize_LineSearchRoots.def_readwrite("success", &Lielab::optimize::LineSearchRoots::success);
    // Lielab_optimize_LineSearchRoots.def_readwrite("message", &Lielab::optimize::LineSearchRoots::message);
    // Lielab_optimize_LineSearchRoots.def(py::init());
    // Lielab_optimize_LineSearchRoots.def("__call__", &Lielab::optimize::LineSearchRoots::operator());
}
