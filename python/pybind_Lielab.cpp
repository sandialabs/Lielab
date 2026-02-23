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

#include "lielab/pybind_domain.hpp"
#include "lielab/pybind_functions.hpp"
#include "lielab/pybind_integrate.hpp"
#include "lielab/pybind_optimize.hpp"
#include "lielab/pybind_testing.hpp"
#include "lielab/pybind_utils.hpp"

namespace py = pybind11;

const bool with_operators = false;

PYBIND11_MODULE(cppLielab, m)
{
    m.doc() = "Lielab Python plugin";
    m.attr("__author__") = Lielab::AUTHOR;
    m.attr("__contact__") = Lielab::CONTACT;
    m.attr("__location__") = Lielab::LOCATION;
    m.attr("__version__") = Lielab::VERSION;

    // Root
    m.def("get_simd_info", &Lielab::get_simd_info);
    m.def("get_eigen_info", &Lielab::get_eigen_info);

    m.def("get_pybind11_info", []()
    {
        return "Pybind11 version: " + std::to_string(PYBIND11_VERSION_MAJOR) + "." + std::to_string(PYBIND11_VERSION_MINOR) + "." + std::to_string(PYBIND11_VERSION_PATCH);
    });

    py::module m_domain = m.def_submodule("domain", "The domain submodule.");
    bind_domain(m_domain);

    py::module m_functions = m.def_submodule("functions", "The functions submodule.");
    bind_functions(m_functions);

    py::module m_integrate = m.def_submodule("integrate", "The integrate submodule.");
    bind_integrate(m_integrate);

    py::module m_optimize = m.def_submodule("optimize", "The optimize submodule.");
    bind_optimize(m_optimize);

    py::module m_testing = m.def_submodule("testing", "The testing submodule.");
    bind_testing(m_testing);

    py::module m_utils = m.def_submodule("utils", "The utils submodule.");
    bind_utils(m_utils);
}
