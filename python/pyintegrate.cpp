#include <Lielab.hpp>

#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "pyintegrate.hpp"

#include <string>
#include <sstream>

namespace py = pybind11;

std::string matstr(const Eigen::VectorXd &mat);
std::string matstr(const Eigen::VectorXcd &mat);
std::string matstr(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mat);
std::string matstr(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> &mat);
std::string matstr(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &mat);
std::string matstr(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &mat);

void bind_integrate(py::module &m_integrate)
{
    m_integrate.def("get_butcher_tableau", &Lielab::integrate::get_butcher_tableau);

    auto Lielab_integrate_Coefficients = py::enum_<Lielab::integrate::Coefficients>(m_integrate, "Coefficients");
    Lielab_integrate_Coefficients.value("E1", Lielab::integrate::Coefficients::E1);
    Lielab_integrate_Coefficients.value("RK3", Lielab::integrate::Coefficients::RK3);
    Lielab_integrate_Coefficients.value("RK4a", Lielab::integrate::Coefficients::RK4a);
    Lielab_integrate_Coefficients.value("RK4b", Lielab::integrate::Coefficients::RK4b);
    Lielab_integrate_Coefficients.value("RK5a", Lielab::integrate::Coefficients::RK5a);
    Lielab_integrate_Coefficients.value("RK5b", Lielab::integrate::Coefficients::RK5b);
    Lielab_integrate_Coefficients.value("RKF12a", Lielab::integrate::Coefficients::RKF12a);
    Lielab_integrate_Coefficients.value("RKF12b", Lielab::integrate::Coefficients::RKF12b);
    Lielab_integrate_Coefficients.value("RKF23a", Lielab::integrate::Coefficients::RKF23a);
    Lielab_integrate_Coefficients.value("RKF23b", Lielab::integrate::Coefficients::RKF23b);
    Lielab_integrate_Coefficients.value("RKF34a", Lielab::integrate::Coefficients::RKF34a);
    Lielab_integrate_Coefficients.value("RKF34b", Lielab::integrate::Coefficients::RKF34b);
    Lielab_integrate_Coefficients.value("RKF45a", Lielab::integrate::Coefficients::RKF45a);
    Lielab_integrate_Coefficients.value("RKF45b", Lielab::integrate::Coefficients::RKF45b);
    Lielab_integrate_Coefficients.value("RKF56", Lielab::integrate::Coefficients::RKF56);
    Lielab_integrate_Coefficients.value("RKF67", Lielab::integrate::Coefficients::RKF67);
    Lielab_integrate_Coefficients.value("RKF78", Lielab::integrate::Coefficients::RKF78);
    Lielab_integrate_Coefficients.value("RKF8", Lielab::integrate::Coefficients::RKF8);
    Lielab_integrate_Coefficients.value("RKDP54_7M", Lielab::integrate::Coefficients::RKDP54_7M);
    Lielab_integrate_Coefficients.value("RKV65e", Lielab::integrate::Coefficients::RKV65e);
    Lielab_integrate_Coefficients.value("RKV65r", Lielab::integrate::Coefficients::RKV65r);
    Lielab_integrate_Coefficients.value("RKV76e", Lielab::integrate::Coefficients::RKV76e);
    Lielab_integrate_Coefficients.value("RKV76r", Lielab::integrate::Coefficients::RKV76r);
    Lielab_integrate_Coefficients.value("RKV87e", Lielab::integrate::Coefficients::RKV87e);
    Lielab_integrate_Coefficients.value("RKV87r", Lielab::integrate::Coefficients::RKV87r);
    Lielab_integrate_Coefficients.value("RKV98e", Lielab::integrate::Coefficients::RKV98e);
    Lielab_integrate_Coefficients.value("RKV98r", Lielab::integrate::Coefficients::RKV98r);
    Lielab_integrate_Coefficients.value("CG4a", Lielab::integrate::Coefficients::CG4a);
    Lielab_integrate_Coefficients.value("CG5a", Lielab::integrate::Coefficients::CG5a);
    
    auto Lielab_integrate_ODESolution = py::class_<Lielab::integrate::ODESolution>(m_integrate, "ODESolution");
    Lielab_integrate_ODESolution.def_readwrite("success", &Lielab::integrate::ODESolution::success);
    Lielab_integrate_ODESolution.def_readwrite("status", &Lielab::integrate::ODESolution::status);
    Lielab_integrate_ODESolution.def_readwrite("message", &Lielab::integrate::ODESolution::message);
    Lielab_integrate_ODESolution.def_readwrite("time_to_solution", &Lielab::integrate::ODESolution::time_to_solution);
    Lielab_integrate_ODESolution.def_readwrite("residuals", &Lielab::integrate::ODESolution::residuals);
    Lielab_integrate_ODESolution.def_readwrite("debug", &Lielab::integrate::ODESolution::debug);
    Lielab_integrate_ODESolution.def_readwrite("t", &Lielab::integrate::ODESolution::t);
    Lielab_integrate_ODESolution.def_readwrite("y", &Lielab::integrate::ODESolution::y);
    Lielab_integrate_ODESolution.def_readwrite("ybar", &Lielab::integrate::ODESolution::ybar);
    Lielab_integrate_ODESolution.def_readwrite("theta", &Lielab::integrate::ODESolution::theta);
    Lielab_integrate_ODESolution.def_readwrite("thetabar", &Lielab::integrate::ODESolution::thetabar);
    Lielab_integrate_ODESolution.def_readwrite("x0", &Lielab::integrate::ODESolution::x0);
    Lielab_integrate_ODESolution.def_readwrite("p", &Lielab::integrate::ODESolution::p);
    Lielab_integrate_ODESolution.def_readwrite("chunk_size", &Lielab::integrate::ODESolution::chunk_size);
    Lielab_integrate_ODESolution.def_readwrite("current_index", &Lielab::integrate::ODESolution::current_index);
    Lielab_integrate_ODESolution.def(py::init<>());
    Lielab_integrate_ODESolution.def(py::init<const Lielab::integrate::ODESolution&>());
    Lielab_integrate_ODESolution.def(py::init<const size_t>());
    // operator =
    Lielab_integrate_ODESolution.def("copy", &Lielab::integrate::ODESolution::copy);
    Lielab_integrate_ODESolution.def("trim_chunk", &Lielab::integrate::ODESolution::trim_chunk);
    Lielab_integrate_ODESolution.def("add_chunk", &Lielab::integrate::ODESolution::add_chunk);
    Lielab_integrate_ODESolution.def("add_data", py::overload_cast<const double, const Eigen::VectorXd&>(&Lielab::integrate::ODESolution::add_data));
    Lielab_integrate_ODESolution.def("__repr__", [](const Lielab::integrate::ODESolution& self)
        {
            return "<lielab.integrate.ODESolution>";
        });
    Lielab_integrate_ODESolution.def("__str__", [](const Lielab::integrate::ODESolution& self)
        {
            std::stringstream out;
            const std::string success = (self.success) ? "True" : "False";

            out << "  success: " << success << "\n";
            out << "   status: " << std::to_string(self.status) << ": " + self.message << "\n";
            out << "        t: NumPy Array (" << std::to_string(self.t.size()) << ",)\n";
            out << "        y: CompositeManifold (" << std::to_string(self.y.size()) << ",)\n";
            out << "       x0: [" << self.x0.transpose() << "]\n";
            out << "        p: [" << self.p.transpose() << "]\n";
            return out.str();
        });

    
    auto Lielab_integrate_RungeKuttaStatus = py::enum_<Lielab::integrate::RungeKuttaStatus>(m_integrate, "RungeKuttaStatus");
    Lielab_integrate_RungeKuttaStatus.value("SUCCESS", Lielab::integrate::RungeKuttaStatus::SUCCESS);
    Lielab_integrate_RungeKuttaStatus.value("DO_STEP0", Lielab::integrate::RungeKuttaStatus::DO_STEP0);
    Lielab_integrate_RungeKuttaStatus.value("DO_STEP1", Lielab::integrate::RungeKuttaStatus::DO_STEP1);
    Lielab_integrate_RungeKuttaStatus.value("ESTIMATE_ERROR", Lielab::integrate::RungeKuttaStatus::ESTIMATE_ERROR);

    auto Lielab_integrate_RungeKutta = py::class_<Lielab::integrate::RungeKutta>(m_integrate, "RungeKutta");
    Lielab_integrate_RungeKutta.def_readwrite("status", &Lielab::integrate::RungeKutta::status);
    Lielab_integrate_RungeKutta.def_readwrite("message", &Lielab::integrate::RungeKutta::message);
    Lielab_integrate_RungeKutta.def_readwrite("error_estimate", &Lielab::integrate::RungeKutta::error_estimate);
    Lielab_integrate_RungeKutta.def_readwrite("can_variable_step", &Lielab::integrate::RungeKutta::can_variable_step);
    Lielab_integrate_RungeKutta.def_readwrite("order", &Lielab::integrate::RungeKutta::order);
    Lielab_integrate_RungeKutta.def_readwrite("A", &Lielab::integrate::RungeKutta::A);
    Lielab_integrate_RungeKutta.def_readwrite("B", &Lielab::integrate::RungeKutta::B);
    Lielab_integrate_RungeKutta.def_readwrite("Bhat", &Lielab::integrate::RungeKutta::Bhat);
    Lielab_integrate_RungeKutta.def_readwrite("C", &Lielab::integrate::RungeKutta::C);
    Lielab_integrate_RungeKutta.def_readwrite("e", &Lielab::integrate::RungeKutta::e);
    Lielab_integrate_RungeKutta.def_readwrite("n", &Lielab::integrate::RungeKutta::n);
    Lielab_integrate_RungeKutta.def_readwrite("reltol", &Lielab::integrate::RungeKutta::reltol);
    Lielab_integrate_RungeKutta.def_readwrite("abstol", &Lielab::integrate::RungeKutta::abstol);
    Lielab_integrate_RungeKutta.def_readwrite("stage", &Lielab::integrate::RungeKutta::stage);
    Lielab_integrate_RungeKutta.def_readwrite("t0", &Lielab::integrate::RungeKutta::t0);
    Lielab_integrate_RungeKutta.def_readwrite("dt", &Lielab::integrate::RungeKutta::dt);
    Lielab_integrate_RungeKutta.def_readwrite("K", &Lielab::integrate::RungeKutta::K);
    Lielab_integrate_RungeKutta.def_readwrite("next_t", &Lielab::integrate::RungeKutta::next_t);
    Lielab_integrate_RungeKutta.def_readwrite("next_theta", &Lielab::integrate::RungeKutta::next_theta);
    Lielab_integrate_RungeKutta.def_readwrite("next_theta2", &Lielab::integrate::RungeKutta::next_theta2);
    
    Lielab_integrate_RungeKutta.def(py::init());
    Lielab_integrate_RungeKutta.def(py::init<Lielab::integrate::Coefficients>());
    Lielab_integrate_RungeKutta.def("init", &Lielab::integrate::RungeKutta::init);
    Lielab_integrate_RungeKutta.def("step_0", &Lielab::integrate::RungeKutta::step_0);
    Lielab_integrate_RungeKutta.def("step_1", &Lielab::integrate::RungeKutta::step_1);
    Lielab_integrate_RungeKutta.def("postprocess", &Lielab::integrate::RungeKutta::postprocess);
    Lielab_integrate_RungeKutta.def("estimate_error", &Lielab::integrate::RungeKutta::estimate_error);
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
    
    auto Lielab_integrate_IVPMethod = py::enum_<Lielab::integrate::IVPMethod>(m_integrate, "IVPMethod");
    Lielab_integrate_IVPMethod.value("Undefined", Lielab::integrate::IVPMethod::Undefined);
    Lielab_integrate_IVPMethod.value("RungeKutta", Lielab::integrate::IVPMethod::RungeKutta);
    Lielab_integrate_IVPMethod.value("MuntheKaas", Lielab::integrate::IVPMethod::MuntheKaas);

    auto Lielab_integrate_IVPOptions = py::class_<Lielab::integrate::IVPOptions>(m_integrate, "IVPOptions");
    Lielab_integrate_IVPOptions.def(py::init<>());
    Lielab_integrate_IVPOptions.def_readwrite("dt_initial", &Lielab::integrate::IVPOptions::dt_initial);
    Lielab_integrate_IVPOptions.def_readwrite("dt_min", &Lielab::integrate::IVPOptions::dt_min);
    Lielab_integrate_IVPOptions.def_readwrite("dt_max", &Lielab::integrate::IVPOptions::dt_max);
    Lielab_integrate_IVPOptions.def_readwrite("reltol", &Lielab::integrate::IVPOptions::reltol);
    Lielab_integrate_IVPOptions.def_readwrite("abstol", &Lielab::integrate::IVPOptions::abstol);
    Lielab_integrate_IVPOptions.def_readwrite("small", &Lielab::integrate::IVPOptions::small);
    Lielab_integrate_IVPOptions.def_readwrite("large", &Lielab::integrate::IVPOptions::large);
    Lielab_integrate_IVPOptions.def_readwrite("pessimist", &Lielab::integrate::IVPOptions::pessimist);
    Lielab_integrate_IVPOptions.def_readwrite("variable_time_step", &Lielab::integrate::IVPOptions::variable_time_step);

    auto Lielab_integrate_RungeKuttaFlowStatus = py::enum_<Lielab::integrate::RungeKuttaFlowStatus>(m_integrate, "RungeKuttaFlowStatus");
    Lielab_integrate_RungeKuttaFlowStatus.value("ERROR_SMALL_DT", Lielab::integrate::RungeKuttaFlowStatus::ERROR_SMALL_DT);
    Lielab_integrate_RungeKuttaFlowStatus.value("ERROR_NEGATIVE_DT", Lielab::integrate::RungeKuttaFlowStatus::ERROR_NEGATIVE_DT);
    Lielab_integrate_RungeKuttaFlowStatus.value("ERROR_SUCCEEDED_BUT_TOL_THO", Lielab::integrate::RungeKuttaFlowStatus::ERROR_SUCCEEDED_BUT_TOL_THO);
    Lielab_integrate_RungeKuttaFlowStatus.value("ERROR_MAX_ITERATIONS", Lielab::integrate::RungeKuttaFlowStatus::ERROR_MAX_ITERATIONS);
    Lielab_integrate_RungeKuttaFlowStatus.value("ERROR_INPUT", Lielab::integrate::RungeKuttaFlowStatus::ERROR_INPUT);
    Lielab_integrate_RungeKuttaFlowStatus.value("SUCCESS", Lielab::integrate::RungeKuttaFlowStatus::SUCCESS);
    Lielab_integrate_RungeKuttaFlowStatus.value("DO_STEP0", Lielab::integrate::RungeKuttaFlowStatus::DO_STEP0);
    Lielab_integrate_RungeKuttaFlowStatus.value("DO_STEP1", Lielab::integrate::RungeKuttaFlowStatus::DO_STEP1);

    auto Lielab_integrate_RungeKuttaFlow = py::class_<Lielab::integrate::RungeKuttaFlow>(m_integrate, "RungeKuttaFlow");
    Lielab_integrate_RungeKuttaFlow.def_readwrite("small", &Lielab::integrate::RungeKuttaFlow::small);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("large", &Lielab::integrate::RungeKuttaFlow::large);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("pessimist", &Lielab::integrate::RungeKuttaFlow::pessimist);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("dt_min", &Lielab::integrate::RungeKuttaFlow::dt_min);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("dt_max", &Lielab::integrate::RungeKuttaFlow::dt_max);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("dt", &Lielab::integrate::RungeKuttaFlow::dt);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("variable_time_step", &Lielab::integrate::RungeKuttaFlow::variable_time_step);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("has_event", &Lielab::integrate::RungeKuttaFlow::has_event);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("event_current", &Lielab::integrate::RungeKuttaFlow::event_current);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("event_next", &Lielab::integrate::RungeKuttaFlow::event_next);

    Lielab_integrate_RungeKuttaFlow.def(py::init());
    Lielab_integrate_RungeKuttaFlow.def("init", &Lielab::integrate::RungeKuttaFlow::init);
    Lielab_integrate_RungeKuttaFlow.def("step0", &Lielab::integrate::RungeKuttaFlow::step0);
    Lielab_integrate_RungeKuttaFlow.def("step1", &Lielab::integrate::RungeKuttaFlow::step1);
    Lielab_integrate_RungeKuttaFlow.def("stepE", &Lielab::integrate::RungeKuttaFlow::stepE);
    Lielab_integrate_RungeKuttaFlow.def("postprocess", &Lielab::integrate::RungeKuttaFlow::postprocess);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("_tcurrent", &Lielab::integrate::RungeKuttaFlow::_tcurrent);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("_ycurrent", &Lielab::integrate::RungeKuttaFlow::_ycurrent);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("_ynext", &Lielab::integrate::RungeKuttaFlow::_ynext);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("iterations", &Lielab::integrate::RungeKuttaFlow::iterations);
    
    Lielab_integrate_RungeKuttaFlow.def_readwrite("_out", &Lielab::integrate::RungeKuttaFlow::_out);
    Lielab_integrate_RungeKuttaFlow.def_readwrite("method", &Lielab::integrate::RungeKuttaFlow::method);
    Lielab_integrate_RungeKuttaFlow.def("__repr__",
        [](const Lielab::integrate::RungeKuttaFlow & self)
        {
            return "<lielab.integrate.RungeKuttaFlow>";
        });
    Lielab_integrate_RungeKuttaFlow.def("__str__",
        [](const Lielab::integrate::RungeKuttaFlow & self)
        {
            return "<lielab.integrate.RungeKuttaFlow>";
        });
}
