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

#include "pydomain.hpp"
#include "pyfunctions.hpp"
#include "pyintegrate.hpp"
#include "pyutils.hpp"

namespace py = pybind11;

std::string matstr(const Eigen::VectorXd &mat)
{
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

std::string matstr(const Eigen::VectorXcd &mat)
{
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

std::string matstr(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mat)
{
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

std::string matstr(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> &mat)
{
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

std::string matstr(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &mat)
{
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

std::string matstr(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &mat)
{
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

const bool with_operators = false;

PYBIND11_MODULE(cppLielab, m)
{
    m.doc() = "Lielab Python plugin";
    m.attr("__author__") = Lielab::AUTHOR;
    m.attr("__contact__") = Lielab::CONTACT;
    m.attr("__location__") = Lielab::LOCATION;
    m.attr("__version__") = Lielab::VERSION;

    py::module m_domain = m.def_submodule("domain", "The domain submodule.");
    bind_domain(m_domain);

    py::module m_functions = m.def_submodule("functions", "The functions submodule.");
    bind_functions(m_functions);

    py::module m_integrate = m.def_submodule("integrate", "The integrate submodule.");
    bind_integrate(m_integrate);

    py::module m_utils = m.def_submodule("utils", "The utils submodule.");
    bind_utils(m_utils);

}
