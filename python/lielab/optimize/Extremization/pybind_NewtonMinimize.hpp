#ifndef PYLIELAB_OPTIMIZE_ROOTSOLVING_NEWTONMINIMIZE_HPP
#define PYLIELAB_OPTIMIZE_ROOTSOLVING_NEWTONMINIMIZE_HPP

#include <pybind11/pybind11.h>
namespace py = pybind11;

void bind_NewtonMinimize(py::module& m_optimize);

#endif
