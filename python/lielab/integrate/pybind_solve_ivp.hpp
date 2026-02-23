#ifndef PYLIELAB_INTEGRATE_SOLVE_IVP_HPP
#define PYLIELAB_INTEGRATE_SOLVE_IVP_HPP

#include <pybind11/pybind11.h>
namespace py = pybind11;

void bind_solve_ivp(py::module& m_integrate);

#endif
