#ifndef PYLIELAB_INTEGRATE_IVPMETHODS_HPP
#define PYLIELAB_INTEGRATE_IVPMETHODS_HPP

#include <pybind11/pybind11.h>
namespace py = pybind11;

void bind_IVPMethods(py::module& m_integrate);

#endif
