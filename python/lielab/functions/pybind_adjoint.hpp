#ifndef PYLIELAB_FUNCTIONS_ADJOINT_HPP
#define PYLIELAB_FUNCTIONS_ADJOINT_HPP

#include <pybind11/pybind11.h>
namespace py = pybind11;

void bind_adjoint(py::module &m_functions);

#endif