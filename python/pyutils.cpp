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

#include "pyutils.hpp"

namespace py = pybind11;

std::string matstr(const Eigen::VectorXd& mat);
std::string matstr(const Eigen::VectorXcd& mat);
std::string matstr(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mat);
std::string matstr(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& mat);
std::string matstr(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& mat);
std::string matstr(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& mat);

void bind_utils(py::module &m_utils)
{
    m_utils.def("factorial", &Lielab::utils::factorial, "The factorial function");

    auto Lielab_utils_golden = py::class_<Lielab::utils::opt_golden>(m_utils, "opt_golden");
    Lielab_utils_golden.def(py::init());
    Lielab_utils_golden.def("init", &Lielab::utils::opt_golden::init);
    Lielab_utils_golden.def("step", &Lielab::utils::opt_golden::step);
    Lielab_utils_golden.def_readwrite("iterations", &Lielab::utils::opt_golden::iterations);
    Lielab_utils_golden.def_readwrite("max_iterations", &Lielab::utils::opt_golden::max_iterations);
    Lielab_utils_golden.def_readwrite("num_objective_evals", &Lielab::utils::opt_golden::num_objective_evals);
    Lielab_utils_golden.def_readwrite("num_jacobian_evals", &Lielab::utils::opt_golden::num_jacobian_evals);
    Lielab_utils_golden.def_readwrite("num_hessian_evals", &Lielab::utils::opt_golden::num_hessian_evals);
    Lielab_utils_golden.def_readwrite("success", &Lielab::utils::opt_golden::success);
    Lielab_utils_golden.def_readwrite("tolerance", &Lielab::utils::opt_golden::tolerance);
    Lielab_utils_golden.def_readwrite("val_objective", &Lielab::utils::opt_golden::val_objective);
    Lielab_utils_golden.def_readwrite("val_jacobian", &Lielab::utils::opt_golden::val_jacobian);
    Lielab_utils_golden.def_readwrite("val_hessian", &Lielab::utils::opt_golden::val_hessian);
    Lielab_utils_golden.def_readwrite("x0", &Lielab::utils::opt_golden::x0);
    Lielab_utils_golden.def_readwrite("lower", &Lielab::utils::opt_golden::lower);
    Lielab_utils_golden.def_readwrite("upper", &Lielab::utils::opt_golden::upper);
    Lielab_utils_golden.def_readwrite("tau", &Lielab::utils::opt_golden::tau);
    Lielab_utils_golden.def_readwrite("_f1", &Lielab::utils::opt_golden::_f1);
    Lielab_utils_golden.def_readwrite("_f2", &Lielab::utils::opt_golden::_f2);
    Lielab_utils_golden.def_readwrite("_X", &Lielab::utils::opt_golden::_X);
    Lielab_utils_golden.def_readwrite("_X1", &Lielab::utils::opt_golden::_X1);
    Lielab_utils_golden.def_readwrite("_X2", &Lielab::utils::opt_golden::_X2);
    Lielab_utils_golden.def_readwrite("_A", &Lielab::utils::opt_golden::_A);
    Lielab_utils_golden.def_readwrite("_B", &Lielab::utils::opt_golden::_B);
    
    auto Lielab_utils_newton = py::class_<Lielab::utils::newton>(m_utils, "newton");
    Lielab_utils_newton.def(py::init());
    Lielab_utils_newton.def("init", &Lielab::utils::newton::init);
    Lielab_utils_newton.def("step0", &Lielab::utils::newton::step0);
    Lielab_utils_newton.def("step1", &Lielab::utils::newton::step1);
    Lielab_utils_newton.def_readwrite("lower", &Lielab::utils::newton::lower);
    Lielab_utils_newton.def_readwrite("upper", &Lielab::utils::newton::upper);
    Lielab_utils_newton.def_readwrite("f", &Lielab::utils::newton::f);
    Lielab_utils_newton.def_readwrite("fnext", &Lielab::utils::newton::fnext);
    Lielab_utils_newton.def_readwrite("m", &Lielab::utils::newton::m);
    Lielab_utils_newton.def_readwrite("mnext", &Lielab::utils::newton::mnext);

    m_utils.def("concatenate", py::overload_cast<const std::vector<Eigen::VectorXi>&>(&Lielab::utils::concatenate<int>));
    m_utils.def("concatenate", py::overload_cast<const std::vector<Eigen::VectorXd>&>(&Lielab::utils::concatenate<double>));

    m_utils.def("arange", py::overload_cast<const int, const int>(&Lielab::utils::arange<int>));
    m_utils.def("arange", py::overload_cast<const double, const double>(&Lielab::utils::arange<double>));
    m_utils.def("arange", py::overload_cast<const int>(&Lielab::utils::arange<int>));
    m_utils.def("arange", py::overload_cast<const double>(&Lielab::utils::arange<double>));

    m_utils.def("repeat", py::overload_cast<const Eigen::VectorXi&, const ptrdiff_t>(&Lielab::utils::repeat<int>));
    m_utils.def("repeat", py::overload_cast<const Eigen::VectorXd&, const ptrdiff_t>(&Lielab::utils::repeat<double>));

    m_utils.def("tile", py::overload_cast<const Eigen::VectorXi&, const ptrdiff_t>(&Lielab::utils::tile<int>));
    m_utils.def("tile", py::overload_cast<const Eigen::VectorXd&, const ptrdiff_t>(&Lielab::utils::tile<double>));

    m_utils.def("linspace", py::overload_cast<const double, const double, const ptrdiff_t>(&Lielab::utils::linspace<double>));

    m_utils.def("column_stack", py::overload_cast<const std::vector<Eigen::VectorXi>&>(&Lielab::utils::column_stack<int>));
    m_utils.def("column_stack", py::overload_cast<const std::vector<Eigen::VectorXd>&>(&Lielab::utils::column_stack<double>));
}
