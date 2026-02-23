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

#include "pybind_RootSolvingCommon.hpp"

namespace py = pybind11;

void bind_RootSolvingCommon(py::module& m_optimize)
{
    auto Lielab_optimize_RootSolvingMethod = py::native_enum<Lielab::optimize::RootSolvingMethod>(m_optimize, "RootSolvingMethod", "enum.Enum");
    Lielab_optimize_RootSolvingMethod.value("Undefined", Lielab::optimize::RootSolvingMethod::Undefined);
    Lielab_optimize_RootSolvingMethod.value("Newton", Lielab::optimize::RootSolvingMethod::Newton);
    Lielab_optimize_RootSolvingMethod.value("LineSearch", Lielab::optimize::RootSolvingMethod::LineSearch);
    Lielab_optimize_RootSolvingMethod.value("Hybrid", Lielab::optimize::RootSolvingMethod::Hybrid);
    Lielab_optimize_RootSolvingMethod.finalize();

    auto Lielab_optimize_RootSolvingOptions = py::class_<Lielab::optimize::RootSolvingOptions>(m_optimize, "RootSolvingOptions");
    Lielab_optimize_RootSolvingOptions.def(py::init<>());
    Lielab_optimize_RootSolvingOptions.def_readwrite("method", &Lielab::optimize::RootSolvingOptions::method);
    Lielab_optimize_RootSolvingOptions.def_readwrite("tol", &Lielab::optimize::RootSolvingOptions::tol);
    Lielab_optimize_RootSolvingOptions.def_readwrite("dx", &Lielab::optimize::RootSolvingOptions::dx);
    Lielab_optimize_RootSolvingOptions.def_readwrite("max_iterations", &Lielab::optimize::RootSolvingOptions::max_iterations);
    Lielab_optimize_RootSolvingOptions.def_readwrite("contraction_factor", &Lielab::optimize::RootSolvingOptions::contraction_factor);
    Lielab_optimize_RootSolvingOptions.def_readwrite("initial_alpha", &Lielab::optimize::RootSolvingOptions::initial_alpha);
    Lielab_optimize_RootSolvingOptions.def_readwrite("initial_step_size", &Lielab::optimize::RootSolvingOptions::initial_step_size);
    Lielab_optimize_RootSolvingOptions.def_readwrite("sufficient_decrease", &Lielab::optimize::RootSolvingOptions::sufficient_decrease);

    auto Lielab_optimize_EuclideanRootSystem = py::class_<Lielab::optimize::EuclideanRootSystem>(m_optimize, "EuclideanRootSystem");
    Lielab_optimize_EuclideanRootSystem.def_readwrite("lower_bound", &Lielab::optimize::EuclideanRootSystem::lower_bound);
    Lielab_optimize_EuclideanRootSystem.def_readwrite("upper_bound", &Lielab::optimize::EuclideanRootSystem::upper_bound);
    Lielab_optimize_EuclideanRootSystem.def(py::init<Lielab::optimize::EuclideanRootSystem_objective_t>());
    Lielab_optimize_EuclideanRootSystem.def_readwrite("objective", &Lielab::optimize::EuclideanRootSystem::objective);
    Lielab_optimize_EuclideanRootSystem.def_readwrite("jacobian", &Lielab::optimize::EuclideanRootSystem::jacobian);

    m_optimize.def("wrap_with_finite_difference", &Lielab::optimize::wrap_with_finite_difference, py::arg("system"), py::arg("x"), py::arg("dx") = 1.0e-6);
}
