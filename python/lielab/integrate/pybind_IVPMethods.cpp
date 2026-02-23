#include <Lielab.hpp>

#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/native_enum.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "pybind_IVPMethods.hpp"
#include "IVPMethods/pybind_Coefficients.hpp"
#include "IVPMethods/pybind_IVPCommon.hpp"
#include "IVPMethods/pybind_RungeKutta.hpp"
#include "IVPMethods/pybind_CrouchGrossman.hpp"
#include "IVPMethods/pybind_MuntheKaas.hpp"
#include "pybind_solve_ivp.hpp"

#include <string>
#include <sstream>

namespace py = pybind11;

void bind_IVPMethods(py::module& m_integrate)
{
    bind_Coefficients(m_integrate);
    bind_IVPCommon(m_integrate);
    bind_RungeKutta(m_integrate);
    bind_CrouchGrossman(m_integrate);
    bind_MuntheKaas(m_integrate);
    
    bind_solve_ivp(m_integrate);
}
