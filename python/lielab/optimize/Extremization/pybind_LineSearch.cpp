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

#include "pybind_LineSearch.hpp"

namespace py = pybind11;

void bind_LineSearch(py::module& m_optimize)
{
    auto Lielab_optimize_LineSearch = py::class_<Lielab::optimize::LineSearch>(m_optimize, "LineSearch");
    Lielab_optimize_LineSearch.def(py::init());
    Lielab_optimize_LineSearch.def_readwrite("iteration", &Lielab::optimize::LineSearch::iteration);
    Lielab_optimize_LineSearch.def_readwrite("success", &Lielab::optimize::LineSearch::success);
    Lielab_optimize_LineSearch.def_readwrite("message", &Lielab::optimize::LineSearch::message);
    Lielab_optimize_LineSearch.def_readwrite("alpha", &Lielab::optimize::LineSearch::alpha);
    Lielab_optimize_LineSearch.def("__call__", py::overload_cast<const Lielab::optimize::EuclideanExtremizationSystem, const Eigen::VectorXd&, const Eigen::VectorXd&, const double, const Lielab::optimize::ExtremizationOptions>(&Lielab::optimize::LineSearch::operator()));
    Lielab_optimize_LineSearch.def("__repr__", [](const Lielab::optimize::LineSearch& self)
        {
            return "<lielab.optimize.LineSearch>";
        });
    Lielab_optimize_LineSearch.def("__str__", [](const Lielab::optimize::LineSearch& self)
        {
            return "<lielab.optimize.LineSearch>";
        });
}
