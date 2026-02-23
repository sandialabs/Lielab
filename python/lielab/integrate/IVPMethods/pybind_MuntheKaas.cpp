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

#include "pybind_MuntheKaas.hpp"

#include <string>
#include <sstream>

namespace py = pybind11;

void bind_MuntheKaas(py::module& m_integrate)
{
    auto Lielab_integrate_MuntheKaas = py::class_<Lielab::integrate::MuntheKaas>(m_integrate, "MuntheKaas");
    Lielab_integrate_MuntheKaas.def_readwrite("status", &Lielab::integrate::MuntheKaas::status);
    Lielab_integrate_MuntheKaas.def_readwrite("success", &Lielab::integrate::MuntheKaas::success);
    Lielab_integrate_MuntheKaas.def_readwrite("message", &Lielab::integrate::MuntheKaas::message);
    Lielab_integrate_MuntheKaas.def_readwrite("abstol", &Lielab::integrate::MuntheKaas::abstol);
    Lielab_integrate_MuntheKaas.def_readwrite("reltol", &Lielab::integrate::MuntheKaas::reltol);
    Lielab_integrate_MuntheKaas.def_readwrite("error_estimate", &Lielab::integrate::MuntheKaas::error_estimate);
    Lielab_integrate_MuntheKaas.def_readwrite("can_variable_step", &Lielab::integrate::MuntheKaas::can_variable_step);
    Lielab_integrate_MuntheKaas.def_readwrite("implicit", &Lielab::integrate::MuntheKaas::implicit);
    Lielab_integrate_MuntheKaas.def_readwrite("order", &Lielab::integrate::MuntheKaas::order);
    Lielab_integrate_MuntheKaas.def_readwrite("A", &Lielab::integrate::MuntheKaas::A);
    Lielab_integrate_MuntheKaas.def_readwrite("B", &Lielab::integrate::MuntheKaas::B);
    Lielab_integrate_MuntheKaas.def_readwrite("Bhat", &Lielab::integrate::MuntheKaas::Bhat);
    Lielab_integrate_MuntheKaas.def_readwrite("C", &Lielab::integrate::MuntheKaas::C);
    Lielab_integrate_MuntheKaas.def_readwrite("e", &Lielab::integrate::MuntheKaas::e);
    Lielab_integrate_MuntheKaas.def_readwrite("n", &Lielab::integrate::MuntheKaas::n);
    Lielab_integrate_MuntheKaas.def_readwrite("K", &Lielab::integrate::MuntheKaas::K);
    Lielab_integrate_MuntheKaas.def(py::init<Lielab::integrate::RungeKuttaCoefficients>(), py::arg("coefficients") = Lielab::integrate::RungeKuttaCoefficients::RKV87r);
    Lielab_integrate_MuntheKaas.def("estimate_error", &Lielab::integrate::MuntheKaas::estimate_error);
    Lielab_integrate_MuntheKaas.def("__call__", &Lielab::integrate::MuntheKaas::operator());
    Lielab_integrate_MuntheKaas.def("__repr__",
        [](const Lielab::integrate::MuntheKaas& self)
        {
            return "<lielab.integrate.MuntheKaas>";
        });
    Lielab_integrate_MuntheKaas.def("__str__",
        [](const Lielab::integrate::MuntheKaas& self)
        {
            return "<lielab.integrate.MuntheKaas>";
        });

    auto Lielab_integrate_MuntheKaasFlow = py::class_<Lielab::integrate::MuntheKaasFlow>(m_integrate, "MuntheKaasFlow");
    Lielab_integrate_MuntheKaasFlow.def_readwrite("status", &Lielab::integrate::MuntheKaasFlow::status);
    Lielab_integrate_MuntheKaasFlow.def_readwrite("success", &Lielab::integrate::MuntheKaasFlow::success);
    Lielab_integrate_MuntheKaasFlow.def_readwrite("message", &Lielab::integrate::MuntheKaasFlow::message);
    Lielab_integrate_MuntheKaasFlow.def_readwrite("tolerance_not_met", &Lielab::integrate::MuntheKaasFlow::tolerance_not_met);
    Lielab_integrate_MuntheKaasFlow.def_readwrite("iterations", &Lielab::integrate::MuntheKaasFlow::iterations);
    Lielab_integrate_MuntheKaasFlow.def(py::init());
    Lielab_integrate_MuntheKaasFlow.def_readwrite("_ynext", &Lielab::integrate::MuntheKaasFlow::_ynext);
    Lielab_integrate_MuntheKaasFlow.def("__call__", &Lielab::integrate::MuntheKaasFlow::operator());
    Lielab_integrate_MuntheKaasFlow.def("__repr__",
        [](const Lielab::integrate::MuntheKaasFlow& self)
        {
            return "<lielab.integrate.MuntheKaasFlow>";
        });
    Lielab_integrate_MuntheKaasFlow.def("__str__",
        [](const Lielab::integrate::MuntheKaasFlow& self)
        {
            return "<lielab.integrate.MuntheKaasFlow>";
        });
}
