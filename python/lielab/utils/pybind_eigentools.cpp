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

#include "pybind_eigentools.hpp"

namespace py = pybind11;

void bind_eigentools(py::module &m_utils)
{
    m_utils.def("concatenate", py::overload_cast<const std::vector<Eigen::VectorXi>&>(&Lielab::utils::concatenate<int>));
    m_utils.def("concatenate", py::overload_cast<const std::vector<Eigen::VectorXd>&>(&Lielab::utils::concatenate<double>));

    m_utils.def("arange", py::overload_cast<const int, const int>(&Lielab::utils::arange<int>));
    m_utils.def("arange", py::overload_cast<const double, const double>(&Lielab::utils::arange<double>));
    m_utils.def("arange", py::overload_cast<const int>(&Lielab::utils::arange<int>));
    m_utils.def("arange", py::overload_cast<const double>(&Lielab::utils::arange<double>));

    m_utils.def("repeat", py::overload_cast<const Eigen::VectorXi&, const int>(&Lielab::utils::repeat<int>));
    m_utils.def("repeat", py::overload_cast<const Eigen::VectorXd&, const int>(&Lielab::utils::repeat<double>));

    m_utils.def("tile", py::overload_cast<const Eigen::VectorXi&, const int>(&Lielab::utils::tile<int>));
    m_utils.def("tile", py::overload_cast<const Eigen::VectorXd&, const int>(&Lielab::utils::tile<double>));

    m_utils.def("linspace", py::overload_cast<const double, const double, const int>(&Lielab::utils::linspace<double>));
    m_utils.def("logspace", py::overload_cast<const double, const double, const int>(&Lielab::utils::logspace<double>));

    m_utils.def("column_stack", py::overload_cast<const std::vector<Eigen::VectorXi>&>(&Lielab::utils::column_stack<int>));
    m_utils.def("column_stack", py::overload_cast<const std::vector<Eigen::VectorXd>&>(&Lielab::utils::column_stack<double>));

    m_utils.def("horizontal_stack", py::overload_cast<const std::vector<Eigen::MatrixXd>&>(&Lielab::utils::horizontal_stack<double>));
    m_utils.def("vertical_stack", py::overload_cast<const std::vector<Eigen::MatrixXd>&>(&Lielab::utils::vertical_stack<double>));

    m_utils.def("linear_interpolate", py::overload_cast<const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::MatrixXd&>(&Lielab::utils::linear_interpolate<double>));
}
