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

#include "pybind_Cayley.hpp"

namespace py = pybind11;

template <typename g>
void bind_cay(py::module& m_functions)
{
    using G = Lielab::domain::LieIII<g>;
    m_functions.def("cay", &Lielab::functions::cay<g>);
    m_functions.def("cayinv", &Lielab::functions::cayinv<G>);
    m_functions.def("dcay", py::overload_cast<const g&, const g&>(&Lielab::functions::dcay<g>));
    m_functions.def("dcayinv", py::overload_cast<const g&, const g&>(&Lielab::functions::dcayinv<g>));

    m_functions.def("cay2", &Lielab::functions::cay2<g>);
}

template <typename... gset>
void do_bind_cay(py::module& m_functions, Lielab::domain::TypeList<gset...>)
{
    (bind_cay<gset>(m_functions), ...);
}

void bind_Cayley(py::module& m_functions)
{
    using namespace Lielab::domain;

    do_bind_cay(m_functions, LieAlgebras{});

    // m_functions.def("dcay", py::overload_cast<const cn&>(&Lielab::functions::dcay<cn>));
    // m_functions.def("dcay", py::overload_cast<const glc&>(&Lielab::functions::dcay<glc>));
    // m_functions.def("dcay", py::overload_cast<const glr&>(&Lielab::functions::dcay<glr>));
    // m_functions.def("dcay", py::overload_cast<const rn&>(&Lielab::functions::dcay<rn>));
    // m_functions.def("dcay", py::overload_cast<const se&>(&Lielab::functions::dcay<se>));
    // m_functions.def("dcay", py::overload_cast<const so&>(&Lielab::functions::dcay<so>));
    // m_functions.def("dcay", py::overload_cast<const sp&>(&Lielab::functions::dcay<sp>));
    // m_functions.def("dcay", py::overload_cast<const su&>(&Lielab::functions::dcay<su>));
    // m_functions.def("dcay", py::overload_cast<const CompositeAlgebra&>(&Lielab::functions::dcay<CompositeAlgebra>));

    // m_functions.def("dcayinv", py::overload_cast<const cn&>(&Lielab::functions::dcayinv<cn>));
    // m_functions.def("dcayinv", py::overload_cast<const glc&>(&Lielab::functions::dcayinv<glc>));
    // m_functions.def("dcayinv", py::overload_cast<const glr&>(&Lielab::functions::dcayinv<glr>));
    // m_functions.def("dcayinv", py::overload_cast<const rn&>(&Lielab::functions::dcayinv<rn>));
    // m_functions.def("dcayinv", py::overload_cast<const se&>(&Lielab::functions::dcayinv<se>));
    // m_functions.def("dcayinv", py::overload_cast<const so&>(&Lielab::functions::dcayinv<so>));
    // m_functions.def("dcayinv", py::overload_cast<const sp&>(&Lielab::functions::dcayinv<sp>));
    // m_functions.def("dcayinv", py::overload_cast<const su&>(&Lielab::functions::dcayinv<su>));
    // m_functions.def("dcayinv", py::overload_cast<const CompositeAlgebra&>(&Lielab::functions::dcayinv<CompositeAlgebra>));
}
