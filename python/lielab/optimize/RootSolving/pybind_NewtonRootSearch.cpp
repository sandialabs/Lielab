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

#include "pybind_NewtonRootSearch.hpp"

namespace py = pybind11;

void bind_NewtonRootSearch(py::module& m_optimize)
{
    auto Lielab_optimize_NewtonRootSearch = py::class_<Lielab::optimize::NewtonRootSearch>(m_optimize, "NewtonRootSearch");
    Lielab_optimize_NewtonRootSearch.def_readwrite("iteration", &Lielab::optimize::NewtonRootSearch::iteration);
    Lielab_optimize_NewtonRootSearch.def_readwrite("success", &Lielab::optimize::NewtonRootSearch::success);
    Lielab_optimize_NewtonRootSearch.def_readwrite("message", &Lielab::optimize::NewtonRootSearch::message);
    Lielab_optimize_NewtonRootSearch.def_readwrite("num_objective_evals", &Lielab::optimize::NewtonRootSearch::num_objective_evals);
    Lielab_optimize_NewtonRootSearch.def(py::init());
    Lielab_optimize_NewtonRootSearch.def("__call__", &Lielab::optimize::NewtonRootSearch::operator());
    Lielab_optimize_NewtonRootSearch.def("__repr__", [](const Lielab::optimize::NewtonRootSearch& self)
        {
            return "<lielab.optimize.NewtonRootSearch>";
        });
    Lielab_optimize_NewtonRootSearch.def("__str__", [](const Lielab::optimize::NewtonRootSearch& self)
        {
            return "<lielab.optimize.NewtonRootSearch>";
        });
}
