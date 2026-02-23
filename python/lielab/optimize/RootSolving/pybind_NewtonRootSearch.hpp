#ifndef PYLIELAB_OPTIMIZE_ROOTSOLVING_NEWTONROOTSEARCH_HPP
#define PYLIELAB_OPTIMIZE_ROOTSOLVING_NEWTONROOTSEARCH_HPP

#include <pybind11/pybind11.h>
namespace py = pybind11;

void bind_NewtonRootSearch(py::module& m_optimize);

#endif
