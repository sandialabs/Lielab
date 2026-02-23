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

#include "pybind_ExtremizationCommon.hpp"

namespace py = pybind11;

void bind_ExtremizationCommon(py::module& m_optimize)
{
    auto Lielab_optimize_ExtremizationMethod = py::native_enum<Lielab::optimize::ExtremizationMethod>(m_optimize, "ExtremizationMethod", "enum.Enum");
    Lielab_optimize_ExtremizationMethod.value("Undefined", Lielab::optimize::ExtremizationMethod::Undefined);
    Lielab_optimize_ExtremizationMethod.value("LineSearch", Lielab::optimize::ExtremizationMethod::LineSearch);
    Lielab_optimize_ExtremizationMethod.value("Newton", Lielab::optimize::ExtremizationMethod::Newton);
    Lielab_optimize_ExtremizationMethod.value("Golden", Lielab::optimize::ExtremizationMethod::Golden);
    Lielab_optimize_ExtremizationMethod.finalize();

    auto Lielab_optimize_ExtremizationOptions = py::class_<Lielab::optimize::ExtremizationOptions>(m_optimize, "ExtremizationOptions");
    Lielab_optimize_ExtremizationOptions.def(py::init<>());
    Lielab_optimize_ExtremizationOptions.def_readwrite("method", &Lielab::optimize::ExtremizationOptions::method);
    Lielab_optimize_ExtremizationOptions.def_readwrite("abstol", &Lielab::optimize::ExtremizationOptions::abstol);
    Lielab_optimize_ExtremizationOptions.def_readwrite("reltol", &Lielab::optimize::ExtremizationOptions::reltol);
    Lielab_optimize_ExtremizationOptions.def_readwrite("max_iterations", &Lielab::optimize::ExtremizationOptions::max_iterations);
    Lielab_optimize_ExtremizationOptions.def_readwrite("contraction_factor", &Lielab::optimize::ExtremizationOptions::contraction_factor);
    Lielab_optimize_ExtremizationOptions.def_readwrite("initial_alpha", &Lielab::optimize::ExtremizationOptions::initial_alpha);
    Lielab_optimize_ExtremizationOptions.def_readwrite("initial_step_size", &Lielab::optimize::ExtremizationOptions::initial_step_size);
    Lielab_optimize_ExtremizationOptions.def_readwrite("sufficient_decrease", &Lielab::optimize::ExtremizationOptions::sufficient_decrease);

    auto Lielab_optimize_EuclideanExtremizationSystem = py::class_<Lielab::optimize::EuclideanExtremizationSystem>(m_optimize, "EuclideanExtremizationSystem");
    Lielab_optimize_EuclideanExtremizationSystem.def_readwrite("lower_bound", &Lielab::optimize::EuclideanExtremizationSystem::lower_bound);
    Lielab_optimize_EuclideanExtremizationSystem.def_readwrite("upper_bound", &Lielab::optimize::EuclideanExtremizationSystem::upper_bound);
    Lielab_optimize_EuclideanExtremizationSystem.def(py::init<Lielab::optimize::EuclideanExtremizationSystem_objective_t>());
    Lielab_optimize_EuclideanExtremizationSystem.def_readwrite("objective", &Lielab::optimize::EuclideanExtremizationSystem::objective);
    Lielab_optimize_EuclideanExtremizationSystem.def_readwrite("jacobian", &Lielab::optimize::EuclideanExtremizationSystem::jacobian);
    Lielab_optimize_EuclideanExtremizationSystem.def_readwrite("hessian", &Lielab::optimize::EuclideanExtremizationSystem::hessian);
}
