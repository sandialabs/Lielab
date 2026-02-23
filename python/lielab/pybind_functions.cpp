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

#include "pybind_functions.hpp"

#include "functions/pybind_adjoint.hpp"
#include "functions/pybind_Cayley.hpp"
#include "functions/pybind_exponential.hpp"

namespace py = pybind11;

void bind_functions(py::module& m_functions)
{
    using namespace Lielab::domain;

    m_functions.def("left_Lie_group_action", &Lielab::functions::left_Lie_group_action, "Default action by left product.");
    m_functions.def("right_Lie_group_action", &Lielab::functions::right_Lie_group_action, "Action by right product.");
    // m_functions.def("copair", &Lielab::functions::copair<glr>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<rn>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<so>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<sp>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<su>, "The copair function.");
    
    bind_adjoint(m_functions);
    bind_Cayley(m_functions);
    
    // m_functions.def("Killing", &Lielab::functions::Killing<cn>, "The Killing function.");
    // m_functions.def("Killing", &Lielab::functions::Killing<glr>, "The Killing function.");
    // m_functions.def("Killing", &Lielab::functions::Killing<glc>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<rn>, "The Killing function.");
    // m_functions.def("Killing", &Lielab::functions::Killing<se>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<so>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<sp>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<su>, "The Killing function."); // TODO: Might be wrong
    m_functions.def("Killingform", &Lielab::functions::Killingform<cn>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<glr>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<glc>, "The Killingform function.");
    // m_functions.def("Killingform", &Lielab::functions::Killingform<se>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<rn>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<so>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<sp>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<su>, "The Killingform function."); // TODO: Might be wrong
    
    bind_exponential(m_functions);

    m_functions.def("pair", &Lielab::functions::pair<cn>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<glr>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<glc>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<rn>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<se>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<so>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<sp>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<su>, "The pair function.");
}
