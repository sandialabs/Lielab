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

#include "pybind_HybridRootSearch.hpp"

namespace py = pybind11;

void bind_HybridRootSearch(py::module& m_optimize)
{
    auto Lielab_optimize_HybridRootSearch = py::class_<Lielab::optimize::HybridRootSearch>(m_optimize, "HybridRootSearch");
    Lielab_optimize_HybridRootSearch.def_readwrite("iterations", &Lielab::optimize::HybridRootSearch::iterations);
    Lielab_optimize_HybridRootSearch.def_readwrite("success", &Lielab::optimize::HybridRootSearch::success);
    Lielab_optimize_HybridRootSearch.def_readwrite("message", &Lielab::optimize::HybridRootSearch::message);
    Lielab_optimize_HybridRootSearch.def(py::init());
    Lielab_optimize_HybridRootSearch.def("__call__", &Lielab::optimize::HybridRootSearch::operator());
    Lielab_optimize_HybridRootSearch.def("__repr__", [](const Lielab::optimize::HybridRootSearch& self)
        {
            return "<lielab.optimize.HybridRootSearch>";
        });
    Lielab_optimize_HybridRootSearch.def("__str__", [](const Lielab::optimize::HybridRootSearch& self)
        {
            return "<lielab.optimize.HybridRootSearch>";
        });
}
