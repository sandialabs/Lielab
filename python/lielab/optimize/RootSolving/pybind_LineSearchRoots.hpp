#ifndef PYLIELAB_OPTIMIZE_ROOTSOLVING_LINESEARCHROOTS_HPP
#define PYLIELAB_OPTIMIZE_ROOTSOLVING_LINESEARCHROOTS_HPP

#include <pybind11/pybind11.h>
namespace py = pybind11;

void bind_LineSearchRoots(py::module& m_optimize);

#endif
