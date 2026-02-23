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

#include "pybind_assertions.hpp"

namespace py = pybind11;

void bind_assertions(py::module& m_testing)
{
    using namespace Lielab::domain;

    m_testing.def("check_topology", py::overload_cast<const CompositeManifold&, const CompositeManifold&>(Lielab::testing::check_topology));
    m_testing.def("check_topology", py::overload_cast<const CompositeGroup&, const CompositeGroup&>(Lielab::testing::check_topology));
    m_testing.def("check_topology", py::overload_cast<const CompositeAlgebra&, const CompositeAlgebra&>(Lielab::testing::check_topology));

    m_testing.def("check_almost_equal_tol", py::overload_cast<const double, const double, const double, const double>(Lielab::testing::check_almost_equal_tol), py::arg("a"), py::arg("b"), py::arg("abstol") = 1.0e-14, py::arg("reltol") = 1.0e-14);
    m_testing.def("check_almost_equal_tol", py::overload_cast<const Eigen::MatrixXd&, const Eigen::MatrixXd&, const double, const double>(Lielab::testing::check_almost_equal_tol), py::arg("a"), py::arg("b"), py::arg("abstol") = 1.0e-14, py::arg("reltol") = 1.0e-14);
    m_testing.def("check_almost_equal_tol", py::overload_cast<const Eigen::MatrixXcd&, const Eigen::MatrixXcd&, const double, const double>(Lielab::testing::check_almost_equal_tol), py::arg("a"), py::arg("b"), py::arg("abstol") = 1.0e-14, py::arg("reltol") = 1.0e-14);
    m_testing.def("check_almost_equal_tol", py::overload_cast<const CompositeManifold&, const CompositeManifold&, const double, const double>(Lielab::testing::check_almost_equal_tol), py::arg("a"), py::arg("b"), py::arg("abstol") = 1.0e-14, py::arg("reltol") = 1.0e-14);
    m_testing.def("check_almost_equal_tol", py::overload_cast<const CompositeGroup&, const CompositeGroup&, const double, const double>(Lielab::testing::check_almost_equal_tol), py::arg("a"), py::arg("b"), py::arg("abstol") = 1.0e-14, py::arg("reltol") = 1.0e-14);
    m_testing.def("check_almost_equal_tol", py::overload_cast<const CompositeAlgebra&, const CompositeAlgebra&, const double, const double>(Lielab::testing::check_almost_equal_tol), py::arg("a"), py::arg("b"), py::arg("abstol") = 1.0e-14, py::arg("reltol") = 1.0e-14);

    m_testing.def("check_almost_equal_nulp", py::overload_cast<const double, const double, const int, const bool>(Lielab::testing::check_almost_equal_nulp), py::arg("a"), py::arg("b"), py::arg("nulp") = 1, py::arg("gate") = false);
    m_testing.def("check_almost_equal_nulp", py::overload_cast<const Eigen::MatrixXd&, const Eigen::MatrixXd&, const int, const bool>(Lielab::testing::check_almost_equal_nulp), py::arg("a"), py::arg("b"), py::arg("nulp") = 1, py::arg("gate") = false);
    m_testing.def("check_almost_equal_nulp", py::overload_cast<const Eigen::MatrixXcd&, const Eigen::MatrixXcd&, const int, const bool>(Lielab::testing::check_almost_equal_nulp), py::arg("a"), py::arg("b"), py::arg("nulp") = 1, py::arg("gate") = false);
    m_testing.def("check_almost_equal_nulp", py::overload_cast<const CompositeManifold&, const CompositeManifold&, const int, const bool>(Lielab::testing::check_almost_equal_nulp), py::arg("a"), py::arg("b"), py::arg("nulp") = 1, py::arg("gate") = false);
    m_testing.def("check_almost_equal_nulp", py::overload_cast<const CompositeGroup&, const CompositeGroup&, const int, const bool>(Lielab::testing::check_almost_equal_nulp), py::arg("a"), py::arg("b"), py::arg("nulp") = 1, py::arg("gate") = false);
    m_testing.def("check_almost_equal_nulp", py::overload_cast<const CompositeAlgebra&, const CompositeAlgebra&, const int, const bool>(Lielab::testing::check_almost_equal_nulp), py::arg("a"), py::arg("b"), py::arg("nulp") = 1, py::arg("gate") = false);
}
