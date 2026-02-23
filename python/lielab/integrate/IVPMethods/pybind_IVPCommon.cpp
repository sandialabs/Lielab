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

#include "pybind_IVPCommon.hpp"

#include <string>
#include <sstream>

namespace py = pybind11;

void bind_IVPCommon(py::module& m_integrate)
{
    auto Lielab_integrate_IVPMethod = py::native_enum<Lielab::integrate::IVPMethod>(m_integrate, "IVPMethod", "enum.Enum");
    Lielab_integrate_IVPMethod.value("Undefined", Lielab::integrate::IVPMethod::Undefined);
    Lielab_integrate_IVPMethod.value("RungeKutta", Lielab::integrate::IVPMethod::RungeKutta);
    Lielab_integrate_IVPMethod.value("CrouchGrossman", Lielab::integrate::IVPMethod::CrouchGrossman);
    Lielab_integrate_IVPMethod.value("MuntheKaas", Lielab::integrate::IVPMethod::MuntheKaas);
    Lielab_integrate_IVPMethod.finalize();

    auto Lielab_integrate_IVPOptions = py::class_<Lielab::integrate::IVPOptions>(m_integrate, "IVPOptions");
    Lielab_integrate_IVPOptions.def(py::init<>());
    Lielab_integrate_IVPOptions.def_readwrite("method", &Lielab::integrate::IVPOptions::method);
    Lielab_integrate_IVPOptions.def_readwrite("dt", &Lielab::integrate::IVPOptions::dt);
    Lielab_integrate_IVPOptions.def_readwrite("dt_min", &Lielab::integrate::IVPOptions::dt_min);
    Lielab_integrate_IVPOptions.def_readwrite("dt_max", &Lielab::integrate::IVPOptions::dt_max);
    Lielab_integrate_IVPOptions.def_readwrite("variable_time_step", &Lielab::integrate::IVPOptions::variable_time_step);
    Lielab_integrate_IVPOptions.def_readwrite("reltol", &Lielab::integrate::IVPOptions::reltol);
    Lielab_integrate_IVPOptions.def_readwrite("abstol", &Lielab::integrate::IVPOptions::abstol);
    Lielab_integrate_IVPOptions.def_readwrite("small", &Lielab::integrate::IVPOptions::small);
    Lielab_integrate_IVPOptions.def_readwrite("large", &Lielab::integrate::IVPOptions::large);
    Lielab_integrate_IVPOptions.def_readwrite("pessimist", &Lielab::integrate::IVPOptions::pessimist);
    Lielab_integrate_IVPOptions.def_readwrite("max_iterations", &Lielab::integrate::IVPOptions::max_iterations);
    Lielab_integrate_IVPOptions.def_readwrite("coefficients", &Lielab::integrate::IVPOptions::coefficients);
    Lielab_integrate_IVPOptions.def_readwrite("rebase_every_step", &Lielab::integrate::IVPOptions::rebase_every_step);
    Lielab_integrate_IVPOptions.def_readwrite("crouch_grossman_coefficients", &Lielab::integrate::IVPOptions::crouch_grossman_coefficients);

    auto Lielab_integrate_IVPStatus = py::native_enum<Lielab::integrate::IVPStatus>(m_integrate, "IVPStatus", "enum.Enum");
    Lielab_integrate_IVPStatus.value("ERROR", Lielab::integrate::IVPStatus::ERROR);
    Lielab_integrate_IVPStatus.value("ERROR_MAX_IERATIONS", Lielab::integrate::IVPStatus::ERROR_MAX_ITERATIONS);
    Lielab_integrate_IVPStatus.value("ERROR_INFS_IN_EVENT", Lielab::integrate::IVPStatus::ERROR_INFS_IN_EVENT);
    Lielab_integrate_IVPStatus.value("ERROR_NANS_IN_EVENT", Lielab::integrate::IVPStatus::ERROR_NANS_IN_EVENT);
    Lielab_integrate_IVPStatus.value("ERROR_INFS_IN_VF", Lielab::integrate::IVPStatus::ERROR_INFS_IN_VF);
    Lielab_integrate_IVPStatus.value("ERROR_NANS_IN_VF", Lielab::integrate::IVPStatus::ERROR_NANS_IN_VF);
    Lielab_integrate_IVPStatus.value("RUNNING", Lielab::integrate::IVPStatus::RUNNING);
    Lielab_integrate_IVPStatus.value("SUCCESS", Lielab::integrate::IVPStatus::SUCCESS);
    Lielab_integrate_IVPStatus.value("SUCCESS_EVENT", Lielab::integrate::IVPStatus::SUCCESS_EVENT);
    Lielab_integrate_IVPStatus.value("SUCCESS_BUT_TOL", Lielab::integrate::IVPStatus::SUCCESS_BUT_TOL);
    Lielab_integrate_IVPStatus.value("SUCCESS_EVENT_BUT_TOL", Lielab::integrate::IVPStatus::SUCCESS_EVENT_BUT_TOL);
    Lielab_integrate_IVPStatus.finalize();

    auto Lielab_integrate_IVPSolution = py::class_<Lielab::integrate::IVPSolution>(m_integrate, "IVPSolution");
    Lielab_integrate_IVPSolution.def_readwrite("success", &Lielab::integrate::IVPSolution::success);
    Lielab_integrate_IVPSolution.def_readwrite("status", &Lielab::integrate::IVPSolution::status);
    Lielab_integrate_IVPSolution.def_readwrite("message", &Lielab::integrate::IVPSolution::message);
    Lielab_integrate_IVPSolution.def_readwrite("time_to_solution", &Lielab::integrate::IVPSolution::time_to_solution);
    Lielab_integrate_IVPSolution.def("to_string", &Lielab::integrate::IVPSolution::to_string);
    Lielab_integrate_IVPSolution.def_readwrite("chunk_size", &Lielab::integrate::IVPSolution::chunk_size);
    Lielab_integrate_IVPSolution.def_readwrite("current_index", &Lielab::integrate::IVPSolution::current_index);
    Lielab_integrate_IVPSolution.def_readwrite("t", &Lielab::integrate::IVPSolution::t);
    Lielab_integrate_IVPSolution.def_readwrite("y", &Lielab::integrate::IVPSolution::y);
    Lielab_integrate_IVPSolution.def_readwrite("ybar", &Lielab::integrate::IVPSolution::ybar);
    Lielab_integrate_IVPSolution.def_readwrite("theta", &Lielab::integrate::IVPSolution::theta);
    Lielab_integrate_IVPSolution.def_readwrite("thetabar", &Lielab::integrate::IVPSolution::thetabar);
    Lielab_integrate_IVPSolution.def_readwrite("debug", &Lielab::integrate::IVPSolution::debug);
    Lielab_integrate_IVPSolution.def(py::init<>());
    Lielab_integrate_IVPSolution.def(py::init<const Lielab::integrate::IVPSolution&>());
    Lielab_integrate_IVPSolution.def(py::init<const size_t>());
    // operator =
    Lielab_integrate_IVPSolution.def("trim_chunk", py::overload_cast<>(&Lielab::integrate::IVPSolution::trim_chunk));
    Lielab_integrate_IVPSolution.def("trim_chunk", py::overload_cast<const ptrdiff_t>(&Lielab::integrate::IVPSolution::trim_chunk));
    Lielab_integrate_IVPSolution.def("add_chunk", &Lielab::integrate::IVPSolution::add_chunk);
    Lielab_integrate_IVPSolution.def("add_data", py::overload_cast<const double, const Eigen::VectorXd&>(&Lielab::integrate::IVPSolution::add_data));
    Lielab_integrate_IVPSolution.def(py::pickle([](const Lielab::integrate::IVPSolution& obj)
        {
            // __getstate__
            return py::make_tuple(obj.success, obj.status, obj.message, obj.time_to_solution,
                                  obj.chunk_size, obj.current_index,
                                  obj.t, obj.y, obj.ybar, obj.theta, obj.thetabar,
                                  obj.debug);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 12)
                throw std::runtime_error("IVPSolution: Invalid state.");

            Lielab::integrate::IVPSolution obj;
            obj.success = t[0].cast<bool>();
            obj.status = t[1].cast<Lielab::integrate::IVPStatus>();
            obj.message = t[2].cast<std::string>();
            obj.time_to_solution = t[3].cast<double>();
            obj.chunk_size = t[4].cast<ptrdiff_t>();
            obj.current_index = t[5].cast<ptrdiff_t>();
            obj.t = t[6].cast<Eigen::VectorXd>();
            obj.y = t[7].cast<std::vector<Lielab::domain::CompositeManifold>>();
            obj.ybar = t[8].cast<Eigen::MatrixXd>();
            obj.theta = t[9].cast<std::vector<Lielab::domain::CompositeAlgebra>>();
            obj.thetabar = t[10].cast<Eigen::MatrixXd>();
            obj.debug = t[11].cast<Eigen::MatrixXd>();

            return obj;
        }));
    Lielab_integrate_IVPSolution.def("__repr__", [](const Lielab::integrate::IVPSolution& self)
        {
            return "<lielab.integrate.IVPSolution>";
        });
    Lielab_integrate_IVPSolution.def("__str__", [](const Lielab::integrate::IVPSolution& self)
        {
            return "<lielab.integrate.IVPSolution>";
        });

    auto Lielab_integrate_EuclideanIVPSystem = py::class_<Lielab::integrate::EuclideanIVPSystem>(m_integrate, "EuclideanIVPSystem");
    Lielab_integrate_EuclideanIVPSystem.def(py::init<Lielab::integrate::EuclideanIVP_vectorfield_t>());
    Lielab_integrate_EuclideanIVPSystem.def_readwrite("event", &Lielab::integrate::EuclideanIVPSystem::event);
    Lielab_integrate_EuclideanIVPSystem.def_readwrite("vectorfield", &Lielab::integrate::EuclideanIVPSystem::vectorfield);

    auto Lielab_integrate_HomogeneousIVPSystem = py::class_<Lielab::integrate::HomogeneousIVPSystem>(m_integrate, "HomogeneousIVPSystem");
    Lielab_integrate_HomogeneousIVPSystem.def(py::init<Lielab::integrate::HomogeneousIVP_generator_t>());
    Lielab_integrate_HomogeneousIVPSystem.def_readwrite("action", &Lielab::integrate::HomogeneousIVPSystem::action);
    Lielab_integrate_HomogeneousIVPSystem.def_readwrite("connection", &Lielab::integrate::HomogeneousIVPSystem::connection);
    Lielab_integrate_HomogeneousIVPSystem.def_readwrite("coordinates", &Lielab::integrate::HomogeneousIVPSystem::coordinates);
    Lielab_integrate_HomogeneousIVPSystem.def_readwrite("event", &Lielab::integrate::HomogeneousIVPSystem::event);
    Lielab_integrate_HomogeneousIVPSystem.def_readwrite("generator", &Lielab::integrate::HomogeneousIVPSystem::generator);

    auto Lielab_integrate_EuclideanClassicHamiltonianIVPSystem = py::class_<Lielab::integrate::EuclideanClassicHamiltonianIVPSystem>(m_integrate, "EuclideanClassicHamiltonianIVPSystem");
    Lielab_integrate_EuclideanClassicHamiltonianIVPSystem.def(py::init<Lielab::integrate::EuclideanClassicHamiltonianIVP_dVdq_t, Lielab::integrate::EuclideanClassicHamiltonianIVP_dTdp_t>());
    Lielab_integrate_EuclideanClassicHamiltonianIVPSystem.def_readwrite("dVdq", &Lielab::integrate::EuclideanClassicHamiltonianIVPSystem::dVdq);
    Lielab_integrate_EuclideanClassicHamiltonianIVPSystem.def_readwrite("dTdp", &Lielab::integrate::EuclideanClassicHamiltonianIVPSystem::dTdp);
    Lielab_integrate_EuclideanClassicHamiltonianIVPSystem.def_readwrite("event", &Lielab::integrate::EuclideanClassicHamiltonianIVPSystem::event);
}
