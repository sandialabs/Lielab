#include <Lielab.hpp>

#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/native_enum.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "pybind_solve_ivp.hpp"

#include <string>
#include <sstream>

namespace py = pybind11;

void bind_solve_ivp(py::module& m_integrate)
{
    m_integrate.def("solve_ivp", py::overload_cast<const Lielab::integrate::EuclideanIVPSystem&, const Eigen::VectorXd&, const Eigen::VectorXd&, const Lielab::integrate::IVPOptions>(&Lielab::integrate::solve_ivp), py::arg("dynamics"), py::arg("tspan"), py::arg("y0"), py::arg("options") = Lielab::integrate::IVPOptions());
    m_integrate.def("solve_ivp", py::overload_cast<const Lielab::integrate::HomogeneousIVPSystem&, const Eigen::VectorXd&, const Lielab::domain::CompositeManifold&, const Lielab::integrate::IVPOptions>(&Lielab::integrate::solve_ivp), py::arg("dynamics"), py::arg("tspan"), py::arg("y0"), py::arg("options") = Lielab::integrate::IVPOptions());
}
