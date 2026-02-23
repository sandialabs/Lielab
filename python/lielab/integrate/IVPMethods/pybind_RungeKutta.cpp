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

#include "pybind_RungeKutta.hpp"

#include <string>
#include <sstream>

namespace py = pybind11;

void bind_RungeKutta(py::module& m_integrate)
{
    auto Lielab_integrate_RungeKutta = py::class_<Lielab::integrate::RungeKutta>(m_integrate, "RungeKutta");
    Lielab_integrate_RungeKutta.def_readwrite("status", &Lielab::integrate::RungeKutta::status);
    Lielab_integrate_RungeKutta.def_readwrite("success", &Lielab::integrate::RungeKutta::success);
    Lielab_integrate_RungeKutta.def_readwrite("message", &Lielab::integrate::RungeKutta::message);
    Lielab_integrate_RungeKutta.def_readwrite("abstol", &Lielab::integrate::RungeKutta::abstol);
    Lielab_integrate_RungeKutta.def_readwrite("reltol", &Lielab::integrate::RungeKutta::reltol);
    Lielab_integrate_RungeKutta.def_readwrite("error_estimate", &Lielab::integrate::RungeKutta::error_estimate);
    Lielab_integrate_RungeKutta.def_readwrite("can_variable_step", &Lielab::integrate::RungeKutta::can_variable_step);
    Lielab_integrate_RungeKutta.def_readwrite("implicit", &Lielab::integrate::RungeKutta::implicit);
    Lielab_integrate_RungeKutta.def_readwrite("order", &Lielab::integrate::RungeKutta::order);
    Lielab_integrate_RungeKutta.def_readwrite("A", &Lielab::integrate::RungeKutta::A);
    Lielab_integrate_RungeKutta.def_readwrite("B", &Lielab::integrate::RungeKutta::B);
    Lielab_integrate_RungeKutta.def_readwrite("Bhat", &Lielab::integrate::RungeKutta::Bhat);
    Lielab_integrate_RungeKutta.def_readwrite("C", &Lielab::integrate::RungeKutta::C);
    Lielab_integrate_RungeKutta.def_readwrite("e", &Lielab::integrate::RungeKutta::e);
    Lielab_integrate_RungeKutta.def_readwrite("n", &Lielab::integrate::RungeKutta::n);
    Lielab_integrate_RungeKutta.def_readwrite("K", &Lielab::integrate::RungeKutta::K);
    Lielab_integrate_RungeKutta.def(py::init<Lielab::integrate::RungeKuttaCoefficients>(), py::arg("coefficients") = Lielab::integrate::RungeKuttaCoefficients::RKV87r);
    Lielab_integrate_RungeKutta.def("estimate_error", &Lielab::integrate::RungeKutta::estimate_error);
    Lielab_integrate_RungeKutta.def("__call__", &Lielab::integrate::RungeKutta::operator());
    Lielab_integrate_RungeKutta.def("__repr__",
        [](const Lielab::integrate::RungeKutta& self)
        {
            return "<lielab.integrate.RungeKutta>";
        });
    Lielab_integrate_RungeKutta.def("__str__",
        [](const Lielab::integrate::RungeKutta& self)
        {
            return "<lielab.integrate.RungeKutta>";
        });

    auto Lielab_integrate_RungeKuttaFlow = py::class_<Lielab::integrate::RungeKuttaFlow>(m_integrate, "RungeKuttaFlow");
    Lielab_integrate_RungeKuttaFlow.def_readwrite("status", &Lielab::integrate::RungeKuttaFlow::status);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("success", &Lielab::integrate::RungeKuttaFlow::success);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("message", &Lielab::integrate::RungeKuttaFlow::message);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("tolerance_not_met", &Lielab::integrate::RungeKuttaFlow::tolerance_not_met);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("iterations", &Lielab::integrate::RungeKuttaFlow::iterations);
    Lielab_integrate_RungeKuttaFlow.def(py::init());
    Lielab_integrate_RungeKuttaFlow.def_readwrite("_ynext", &Lielab::integrate::RungeKuttaFlow::_ynext);
    Lielab_integrate_RungeKuttaFlow.def("__call__", &Lielab::integrate::RungeKuttaFlow::operator());
    Lielab_integrate_RungeKuttaFlow.def("__repr__",
        [](const Lielab::integrate::RungeKuttaFlow& self)
        {
            return "<lielab.integrate.RungeKuttaFlow>";
        });
    Lielab_integrate_RungeKuttaFlow.def("__str__",
        [](const Lielab::integrate::RungeKuttaFlow& self)
        {
            return "<lielab.integrate.RungeKuttaFlow>";
        });
}
