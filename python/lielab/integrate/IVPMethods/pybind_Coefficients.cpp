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

#include "pybind_Coefficients.hpp"

#include <string>
#include <sstream>

namespace py = pybind11;

void bind_Coefficients(py::module& m_integrate)
{
    auto Lielab_integrate_RungeKuttaCoefficients = py::native_enum<Lielab::integrate::RungeKuttaCoefficients>(m_integrate, "RungeKuttaCoefficients", "enum.Enum");
    Lielab_integrate_RungeKuttaCoefficients.value("FE1", Lielab::integrate::RungeKuttaCoefficients::FE1);
    Lielab_integrate_RungeKuttaCoefficients.value("RK3", Lielab::integrate::RungeKuttaCoefficients::RK3);
    Lielab_integrate_RungeKuttaCoefficients.value("RK4a", Lielab::integrate::RungeKuttaCoefficients::RK4a);
    Lielab_integrate_RungeKuttaCoefficients.value("RK4b", Lielab::integrate::RungeKuttaCoefficients::RK4b);
    Lielab_integrate_RungeKuttaCoefficients.value("RK5a", Lielab::integrate::RungeKuttaCoefficients::RK5a);
    Lielab_integrate_RungeKuttaCoefficients.value("RK5b", Lielab::integrate::RungeKuttaCoefficients::RK5b);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF12a", Lielab::integrate::RungeKuttaCoefficients::RKF12a);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF12b", Lielab::integrate::RungeKuttaCoefficients::RKF12b);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF23a", Lielab::integrate::RungeKuttaCoefficients::RKF23a);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF23b", Lielab::integrate::RungeKuttaCoefficients::RKF23b);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF34a", Lielab::integrate::RungeKuttaCoefficients::RKF34a);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF34b", Lielab::integrate::RungeKuttaCoefficients::RKF34b);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF45a", Lielab::integrate::RungeKuttaCoefficients::RKF45a);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF45b", Lielab::integrate::RungeKuttaCoefficients::RKF45b);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF56", Lielab::integrate::RungeKuttaCoefficients::RKF56);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF67", Lielab::integrate::RungeKuttaCoefficients::RKF67);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF78", Lielab::integrate::RungeKuttaCoefficients::RKF78);
    Lielab_integrate_RungeKuttaCoefficients.value("RKF8", Lielab::integrate::RungeKuttaCoefficients::RKF8);
    Lielab_integrate_RungeKuttaCoefficients.value("RKDP54_7M", Lielab::integrate::RungeKuttaCoefficients::RKDP54_7M);
    Lielab_integrate_RungeKuttaCoefficients.value("RKV65e", Lielab::integrate::RungeKuttaCoefficients::RKV65e);
    Lielab_integrate_RungeKuttaCoefficients.value("RKV65r", Lielab::integrate::RungeKuttaCoefficients::RKV65r);
    Lielab_integrate_RungeKuttaCoefficients.value("RKV76e", Lielab::integrate::RungeKuttaCoefficients::RKV76e);
    Lielab_integrate_RungeKuttaCoefficients.value("RKV76r", Lielab::integrate::RungeKuttaCoefficients::RKV76r);
    Lielab_integrate_RungeKuttaCoefficients.value("RKV87e", Lielab::integrate::RungeKuttaCoefficients::RKV87e);
    Lielab_integrate_RungeKuttaCoefficients.value("RKV87r", Lielab::integrate::RungeKuttaCoefficients::RKV87r);
    Lielab_integrate_RungeKuttaCoefficients.value("RKV98e", Lielab::integrate::RungeKuttaCoefficients::RKV98e);
    Lielab_integrate_RungeKuttaCoefficients.value("RKV98r", Lielab::integrate::RungeKuttaCoefficients::RKV98r);
    Lielab_integrate_RungeKuttaCoefficients.value("BE1", Lielab::integrate::RungeKuttaCoefficients::BE1);
    Lielab_integrate_RungeKuttaCoefficients.value("LG2", Lielab::integrate::RungeKuttaCoefficients::LG2);
    Lielab_integrate_RungeKuttaCoefficients.value("LG4", Lielab::integrate::RungeKuttaCoefficients::LG4);
    Lielab_integrate_RungeKuttaCoefficients.value("LG4s", Lielab::integrate::RungeKuttaCoefficients::LG4s);
    Lielab_integrate_RungeKuttaCoefficients.value("LG6", Lielab::integrate::RungeKuttaCoefficients::LG6);
    Lielab_integrate_RungeKuttaCoefficients.value("LG6s", Lielab::integrate::RungeKuttaCoefficients::LG6s);
    Lielab_integrate_RungeKuttaCoefficients.value("Lobatto3A2", Lielab::integrate::RungeKuttaCoefficients::Lobatto3A2);
    Lielab_integrate_RungeKuttaCoefficients.value("Lobatto3A4", Lielab::integrate::RungeKuttaCoefficients::Lobatto3A4);
    Lielab_integrate_RungeKuttaCoefficients.value("Lobatto3A6", Lielab::integrate::RungeKuttaCoefficients::Lobatto3A6);
    Lielab_integrate_RungeKuttaCoefficients.finalize();

    m_integrate.def("get_butcher_tableau", &Lielab::integrate::get_butcher_tableau);

    auto Lielab_integrate_CrouchGrossmanCoefficients = py::native_enum<Lielab::integrate::CrouchGrossmanCoefficients>(m_integrate, "CrouchGrossmanCoefficients", "enum.Enum");
    Lielab_integrate_CrouchGrossmanCoefficients.value("CG23", Lielab::integrate::CrouchGrossmanCoefficients::CG23);
    Lielab_integrate_CrouchGrossmanCoefficients.value("CG4a", Lielab::integrate::CrouchGrossmanCoefficients::CG4a);
    Lielab_integrate_CrouchGrossmanCoefficients.value("CG5a", Lielab::integrate::CrouchGrossmanCoefficients::CG5a);
    Lielab_integrate_CrouchGrossmanCoefficients.finalize();

    m_integrate.def("get_crouch_grossman_coefficients", &Lielab::integrate::get_crouch_grossman_coefficients);
}
