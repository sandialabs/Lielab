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

#include "pybind_CrouchGrossman.hpp"

#include <string>
#include <sstream>

namespace py = pybind11;

void bind_CrouchGrossman(py::module& m_integrate)
{
    auto Lielab_integrate_CrouchGrossman = py::class_<Lielab::integrate::CrouchGrossman>(m_integrate, "CrouchGrossman");
    Lielab_integrate_CrouchGrossman.def_readwrite("status", &Lielab::integrate::CrouchGrossman::status);
    Lielab_integrate_CrouchGrossman.def_readwrite("success", &Lielab::integrate::CrouchGrossman::success);
    Lielab_integrate_CrouchGrossman.def_readwrite("message", &Lielab::integrate::CrouchGrossman::message);
    Lielab_integrate_CrouchGrossman.def_readwrite("abstol", &Lielab::integrate::CrouchGrossman::abstol);
    Lielab_integrate_CrouchGrossman.def_readwrite("reltol", &Lielab::integrate::CrouchGrossman::reltol);
    Lielab_integrate_CrouchGrossman.def_readwrite("error_estimate", &Lielab::integrate::CrouchGrossman::error_estimate);
    Lielab_integrate_CrouchGrossman.def_readwrite("can_variable_step", &Lielab::integrate::CrouchGrossman::can_variable_step);
    Lielab_integrate_CrouchGrossman.def_readwrite("implicit", &Lielab::integrate::CrouchGrossman::implicit);
    Lielab_integrate_CrouchGrossman.def_readwrite("order", &Lielab::integrate::CrouchGrossman::order);
    Lielab_integrate_CrouchGrossman.def_readwrite("A", &Lielab::integrate::CrouchGrossman::A);
    Lielab_integrate_CrouchGrossman.def_readwrite("B", &Lielab::integrate::CrouchGrossman::B);
    Lielab_integrate_CrouchGrossman.def_readwrite("Bhat", &Lielab::integrate::CrouchGrossman::Bhat);
    Lielab_integrate_CrouchGrossman.def_readwrite("C", &Lielab::integrate::CrouchGrossman::C);
    Lielab_integrate_CrouchGrossman.def_readwrite("e", &Lielab::integrate::CrouchGrossman::e);
    Lielab_integrate_CrouchGrossman.def_readwrite("n", &Lielab::integrate::CrouchGrossman::n);
    Lielab_integrate_CrouchGrossman.def_readwrite("K", &Lielab::integrate::CrouchGrossman::K);
    Lielab_integrate_CrouchGrossman.def(py::init<Lielab::integrate::CrouchGrossmanCoefficients>(), py::arg("coefficients") = Lielab::integrate::CrouchGrossmanCoefficients::CG23);
    Lielab_integrate_CrouchGrossman.def("estimate_error", &Lielab::integrate::CrouchGrossman::estimate_error);
    Lielab_integrate_CrouchGrossman.def("__call__", &Lielab::integrate::CrouchGrossman::operator());
    Lielab_integrate_CrouchGrossman.def("__repr__",
        [](const Lielab::integrate::CrouchGrossman& self)
        {
            return "<lielab.integrate.CrouchGrossman>";
        });
    Lielab_integrate_CrouchGrossman.def("__str__",
        [](const Lielab::integrate::CrouchGrossman& self)
        {
            return "<lielab.integrate.CrouchGrossman>";
        });

    auto Lielab_integrate_CrouchGrossmanFlow = py::class_<Lielab::integrate::CrouchGrossmanFlow>(m_integrate, "CrouchGrossmanFlow");
    Lielab_integrate_CrouchGrossmanFlow.def_readwrite("status", &Lielab::integrate::CrouchGrossmanFlow::status);
    Lielab_integrate_CrouchGrossmanFlow.def_readwrite("success", &Lielab::integrate::CrouchGrossmanFlow::success);
    Lielab_integrate_CrouchGrossmanFlow.def_readwrite("message", &Lielab::integrate::CrouchGrossmanFlow::message);
    Lielab_integrate_CrouchGrossmanFlow.def_readwrite("tolerance_not_met", &Lielab::integrate::CrouchGrossmanFlow::tolerance_not_met);
    Lielab_integrate_CrouchGrossmanFlow.def_readwrite("iterations", &Lielab::integrate::CrouchGrossmanFlow::iterations);
    Lielab_integrate_CrouchGrossmanFlow.def(py::init());
    Lielab_integrate_CrouchGrossmanFlow.def_readwrite("_ynext", &Lielab::integrate::CrouchGrossmanFlow::_ynext);
    Lielab_integrate_CrouchGrossmanFlow.def("__call__", &Lielab::integrate::CrouchGrossmanFlow::operator());
    Lielab_integrate_CrouchGrossmanFlow.def("__repr__",
        [](const Lielab::integrate::CrouchGrossmanFlow& self)
        {
            return "<lielab.integrate.CrouchGrossmanFlow>";
        });
    Lielab_integrate_CrouchGrossmanFlow.def("__str__",
        [](const Lielab::integrate::CrouchGrossmanFlow& self)
        {
            return "<lielab.integrate.CrouchGrossmanFlow>";
        });
}
