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

#include "pybind_NewtonMinimize.hpp"

namespace py = pybind11;

void bind_NewtonMinimize(py::module& m_optimize)
{
    auto Lielab_optimize_NewtonMinimize = py::class_<Lielab::optimize::NewtonMinimize>(m_optimize, "NewtonMinimize");
    Lielab_optimize_NewtonMinimize.def_readwrite("iteration", &Lielab::optimize::NewtonMinimize::iteration);
    Lielab_optimize_NewtonMinimize.def_readwrite("success", &Lielab::optimize::NewtonMinimize::success);
    Lielab_optimize_NewtonMinimize.def_readwrite("message", &Lielab::optimize::NewtonMinimize::message);
    Lielab_optimize_NewtonMinimize.def(py::init());
    Lielab_optimize_NewtonMinimize.def("__call__", py::overload_cast<Lielab::optimize::EuclideanExtremizationSystem, const Eigen::VectorXd&, const Lielab::optimize::ExtremizationOptions>(&Lielab::optimize::NewtonMinimize::operator()));
    Lielab_optimize_NewtonMinimize.def("__repr__", [](const Lielab::optimize::NewtonMinimize& self)
        {
            return "<lielab.optimize.NewtonMinimize>";
        });
    Lielab_optimize_NewtonMinimize.def("__str__", [](const Lielab::optimize::NewtonMinimize& self)
        {
            return "<lielab.optimize.NewtonMinimize>";
        });
}
