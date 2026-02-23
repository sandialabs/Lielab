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

#include "pybind_adjoint.hpp"

namespace py = pybind11;

template <typename g>
void bind_ad(py::module& m_functions)
{
    using G = Lielab::domain::LieIII<g>;
    m_functions.def("commutator", &Lielab::functions::commutator<g>);
    m_functions.def("ad_numerical", py::overload_cast<const g&, const g&, const int>(&Lielab::functions::ad_numerical<g>), py::arg("x"), py::arg("y"), py::arg("power") = 1);
    m_functions.def("ad", py::overload_cast<const g&, const g&, const int>(&Lielab::functions::ad<g>), py::arg("x"), py::arg("y"), py::arg("power") = 1);
    m_functions.def("Ad", py::overload_cast<const G&, const g&>(&Lielab::functions::Ad<G>), py::arg("g"), py::arg("x"));

    m_functions.def("coad_numerical", py::overload_cast<const g&, const g&, const int>(&Lielab::functions::coad_numerical<g>), py::arg("x"), py::arg("y"), py::arg("power") = 1);
    m_functions.def("coad", py::overload_cast<const g&, const g&, const int>(&Lielab::functions::coad<g>), py::arg("x"), py::arg("y"), py::arg("power") = 1);
}

template <typename... gset>
void do_bind_ad(py::module& m_functions, Lielab::domain::TypeList<gset...>)
{
    (bind_ad<gset>(m_functions), ...);
}

void bind_adjoint(py::module& m_functions)
{
    using namespace Lielab::domain;

    do_bind_ad(m_functions, LieAlgebras{});

    m_functions.def("ad_numerical", py::overload_cast<const cn&, const int>(&Lielab::functions::ad_numerical<cn>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad_numerical", py::overload_cast<const glc&, const int>(&Lielab::functions::ad_numerical<glc>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad_numerical", py::overload_cast<const glr&, const int>(&Lielab::functions::ad_numerical<glr>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad_numerical", py::overload_cast<const rn&, const int>(&Lielab::functions::ad_numerical<rn>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad_numerical", py::overload_cast<const se&, const int>(&Lielab::functions::ad_numerical<se>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad_numerical", py::overload_cast<const so&, const int>(&Lielab::functions::ad_numerical<so>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad_numerical", py::overload_cast<const sp&, const int>(&Lielab::functions::ad_numerical<sp>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad_numerical", py::overload_cast<const su&, const int>(&Lielab::functions::ad_numerical<su>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad_numerical", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::ad_numerical<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("p") = 1);
    
    m_functions.def("ad", py::overload_cast<const cn&, const int>(&Lielab::functions::ad<cn>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad", py::overload_cast<const glc&, const int>(&Lielab::functions::ad<glc>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad", py::overload_cast<const glr&, const int>(&Lielab::functions::ad<glr>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad", py::overload_cast<const rn&, const int>(&Lielab::functions::ad<rn>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad", py::overload_cast<const se&, const int>(&Lielab::functions::ad<se>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad", py::overload_cast<const so&, const int>(&Lielab::functions::ad<so>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad", py::overload_cast<const sp&, const int>(&Lielab::functions::ad<sp>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad", py::overload_cast<const su&, const int>(&Lielab::functions::ad<su>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("ad", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::ad<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("p") = 1);
    
    m_functions.def("coad_numerical", py::overload_cast<const cn&, const int>(&Lielab::functions::coad_numerical<cn>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad_numerical", py::overload_cast<const glc&, const int>(&Lielab::functions::coad_numerical<glc>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad_numerical", py::overload_cast<const glr&, const int>(&Lielab::functions::coad_numerical<glr>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad_numerical", py::overload_cast<const rn&, const int>(&Lielab::functions::coad_numerical<rn>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad_numerical", py::overload_cast<const se&, const int>(&Lielab::functions::coad_numerical<se>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad_numerical", py::overload_cast<const so&, const int>(&Lielab::functions::coad_numerical<so>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad_numerical", py::overload_cast<const sp&, const int>(&Lielab::functions::coad_numerical<sp>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad_numerical", py::overload_cast<const su&, const int>(&Lielab::functions::coad_numerical<su>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad_numerical", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::coad_numerical<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("p") = 1);
    
    m_functions.def("coad", py::overload_cast<const cn&, const int>(&Lielab::functions::coad<cn>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad", py::overload_cast<const glc&, const int>(&Lielab::functions::coad<glc>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad", py::overload_cast<const glr&, const int>(&Lielab::functions::coad<glr>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad", py::overload_cast<const rn&, const int>(&Lielab::functions::coad<rn>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad", py::overload_cast<const se&, const int>(&Lielab::functions::coad<se>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad", py::overload_cast<const so&, const int>(&Lielab::functions::coad<so>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad", py::overload_cast<const sp&, const int>(&Lielab::functions::coad<sp>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad", py::overload_cast<const su&, const int>(&Lielab::functions::coad<su>), py::arg("a"), py::arg("p") = 1);
    m_functions.def("coad", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::coad<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("p") = 1);

    // m_functions.def("Ad_numerical", [](const cn& a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"));
    // m_functions.def("Ad_numerical", [](const glr& a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"));
    // m_functions.def("Ad_numerical", [](const glc& a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"));
    // m_functions.def("Ad_numerical", [](const rn& a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"));
    // m_functions.def("Ad_numerical", [](const se& a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"));
    // m_functions.def("Ad_numerical", [](const so& a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"));
    // m_functions.def("Ad_numerical", [](const sp& a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"));
    // m_functions.def("Ad_numerical", [](const su& a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"));
    // m_functions.def("Ad", [](const cn& a){return Lielab::functions::Ad(a);}, py::arg("a"));
    // m_functions.def("Ad", [](const glr& a){return Lielab::functions::Ad(a);}, py::arg("a"));
    // m_functions.def("Ad", [](const glc& a){return Lielab::functions::Ad(a);}, py::arg("a"));
    // m_functions.def("Ad", [](const rn& a){return Lielab::functions::Ad(a);}, py::arg("a"));
    // m_functions.def("Ad", [](const se& a){return Lielab::functions::Ad(a);}, py::arg("a"));
    // m_functions.def("Ad", [](const so& a){return Lielab::functions::Ad(a);}, py::arg("a"));
    // m_functions.def("Ad", [](const sp& a){return Lielab::functions::Ad(a);}, py::arg("a"));
    // m_functions.def("Ad", [](const su& a){return Lielab::functions::Ad(a);}, py::arg("a"));
    // m_functions.def("Ad", [](const cn& a, cn& b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"));
    // m_functions.def("Ad", [](const glr& a, glr& b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"));
    // m_functions.def("Ad", [](const glc& a, glc& b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"));
    // m_functions.def("Ad", [](const rn& a, rn& b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"));
    // m_functions.def("Ad", [](const se& a, se& b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"));
    // m_functions.def("Ad", [](const so& a, so& b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"));
    // m_functions.def("Ad", [](const sp& a, sp& b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"));
    // m_functions.def("Ad", [](const su& a, su& b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"));
}
