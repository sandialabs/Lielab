#ifndef PYLIELAB_FUNCTIONS_HPP
#define PYLIELAB_FUNCTIONS_HPP

#include <pybind11/pybind11.h>
namespace py = pybind11;

void bind_functions(py::module &m_domain);

#endif
