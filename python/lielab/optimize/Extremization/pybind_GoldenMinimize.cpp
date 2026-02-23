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

#include "pybind_GoldenMinimize.hpp"

namespace py = pybind11;

void bind_GoldenMinimize(py::module& m_optimize)
{
    auto Lielab_optimize_GoldenMinimize = py::class_<Lielab::optimize::GoldenMinimize>(m_optimize, "GoldenMinimize");
    Lielab_optimize_GoldenMinimize.def(py::init());
    Lielab_optimize_GoldenMinimize.def_readwrite("iteration", &Lielab::optimize::GoldenMinimize::iteration);
    Lielab_optimize_GoldenMinimize.def_readwrite("success", &Lielab::optimize::GoldenMinimize::success);
    Lielab_optimize_GoldenMinimize.def_readwrite("message", &Lielab::optimize::GoldenMinimize::message);
    Lielab_optimize_GoldenMinimize.def("__call__", &Lielab::optimize::GoldenMinimize::operator());
    Lielab_optimize_GoldenMinimize.def("__repr__", [](const Lielab::optimize::GoldenMinimize& self)
        {
            return "<lielab.optimize.GoldenMinimize>";
        });
    Lielab_optimize_GoldenMinimize.def("__str__", [](const Lielab::optimize::GoldenMinimize& self)
        {
            return "<lielab.optimize.GoldenMinimize>";
        });
}
