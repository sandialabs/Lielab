#ifndef PYLIELAB_UTILS_OPTIMIZE_LINESEARCH_HPP
#define PYLIELAB_UTILS_OPTIMIZE_LINESEARCH_HPP

#include <pybind11/pybind11.h>
namespace py = pybind11;

void bind_LineSearch(py::module& m_optimize);

#endif
