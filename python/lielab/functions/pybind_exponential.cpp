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

#include "pybind_exponential.hpp"

namespace py = pybind11;

template <typename g>
void bind_exp(py::module& m_functions)
{
    using Lielab::functions::exp_numerical;
    using Lielab::functions::exp;
    using Lielab::functions::dexp_numerical;
    using Lielab::functions::dexp;
    using Lielab::functions::dexpinv_numerical;
    using Lielab::functions::dexpinv;

    m_functions.def("exp_numerical", &exp_numerical<g>);
    m_functions.def("exp", &exp<g>);
    // m_functions.def("dexp_numerical", py::overload_cast<const g&, const int>(&dexp_numerical<g>), py::arg("x"), py::arg("order") = 5);
    // m_functions.def("dexp", py::overload_cast<const g&, const int>(&dexp<g>), py::arg("x"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", py::overload_cast<const g&, const g&, const int>(&dexp_numerical<g>), py::arg("x"), py::arg("y"), py::arg("order") = 5);
    m_functions.def("dexp", py::overload_cast<const g&, const g&, const int>(&dexp<g>), py::arg("x"), py::arg("y"), py::arg("order") = 5);
    // m_functions.def("dexpinv_numerical", py::overload_cast<const g&, const int>(&dexpinv_numerical<g>), py::arg("x"), py::arg("order") = 5);
    // m_functions.def("dexpinv", py::overload_cast<const g&, const int>(&dexpinv<g>), py::arg("x"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", py::overload_cast<const g&, const g&, const int>(&dexpinv_numerical<g>), py::arg("x"), py::arg("y"), py::arg("order") = 5);
    m_functions.def("dexpinv", py::overload_cast<const g&, const g&, const int>(&dexpinv<g>), py::arg("x"), py::arg("y"), py::arg("order") = 5);
}

template <typename g>
void bind_log(py::module& m_functions)
{
    using Lielab::functions::log_numerical;
    using Lielab::functions::log;
    using Lielab::functions::dlog_numerical;
    using Lielab::functions::dlog;
    using Lielab::functions::dloginv_numerical;
    using Lielab::functions::dloginv;

    using G = Lielab::domain::LieIII<g>;

    m_functions.def("log_numerical", &Lielab::functions::log_numerical<G>);
    m_functions.def("log", &Lielab::functions::log<G>);
    // m_functions.def("dlog_numerical", py::overload_cast<const g&, const int>(&Lielab::functions::dlog_numerical<g>), py::arg("x"), py::arg("order") = 5);
    // m_functions.def("dlog", py::overload_cast<const g&, const int>(&Lielab::functions::dlog<g>), py::arg("x"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", py::overload_cast<const g&, const g&, const int>(&Lielab::functions::dlog_numerical<g>), py::arg("x"), py::arg("y"), py::arg("order") = 5);
    m_functions.def("dlog", py::overload_cast<const g&, const g&, const int>(&Lielab::functions::dlog<g>), py::arg("x"), py::arg("y"), py::arg("order") = 5);
    // m_functions.def("dloginv_numerical", py::overload_cast<const g&, const int>(&Lielab::functions::dloginv_numerical<g>), py::arg("x"), py::arg("order") = 5);
    // m_functions.def("dloginv", py::overload_cast<const g&, const int>(&Lielab::functions::dloginv<g>), py::arg("x"), py::arg("order") = 5);
    m_functions.def("dloginv_numerical", py::overload_cast<const g&, const g&, const int>(&Lielab::functions::dloginv_numerical<g>), py::arg("x"), py::arg("y"), py::arg("order") = 5);
    m_functions.def("dloginv", py::overload_cast<const g&, const g&, const int>(&Lielab::functions::dloginv<g>), py::arg("x"), py::arg("y"), py::arg("order") = 5);
}

template <typename... gset>
void do_bind_exponential(py::module& m_functions, Lielab::domain::TypeList<gset...>)
{
    (bind_exp<gset>(m_functions), ...);
    (bind_log<gset>(m_functions), ...);
}

void bind_exponential(py::module& m_functions)
{
    using namespace Lielab::domain;

    do_bind_exponential(m_functions, LieAlgebras{});

    m_functions.def("dexp_numerical", py::overload_cast<const cn&, const int>(&Lielab::functions::dexp_numerical<cn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", py::overload_cast<const glc&, const int>(&Lielab::functions::dexp_numerical<glc>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", py::overload_cast<const glr&, const int>(&Lielab::functions::dexp_numerical<glr>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", py::overload_cast<const rn&, const int>(&Lielab::functions::dexp_numerical<rn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", py::overload_cast<const se&, const int>(&Lielab::functions::dexp_numerical<se>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", py::overload_cast<const so&, const int>(&Lielab::functions::dexp_numerical<so>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", py::overload_cast<const sp&, const int>(&Lielab::functions::dexp_numerical<sp>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", py::overload_cast<const su&, const int>(&Lielab::functions::dexp_numerical<su>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::dexp_numerical<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("order") = 5);
    
    m_functions.def("dexp", py::overload_cast<const cn&, const int>(&Lielab::functions::dexp<cn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", py::overload_cast<const glc&, const int>(&Lielab::functions::dexp<glc>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", py::overload_cast<const glr&, const int>(&Lielab::functions::dexp<glr>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", py::overload_cast<const rn&, const int>(&Lielab::functions::dexp<rn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", py::overload_cast<const se&, const int>(&Lielab::functions::dexp<se>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", py::overload_cast<const so&, const int>(&Lielab::functions::dexp<so>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", py::overload_cast<const sp&, const int>(&Lielab::functions::dexp<sp>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", py::overload_cast<const su&, const int>(&Lielab::functions::dexp<su>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::dexp<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("order") = 5);
    
    m_functions.def("dexpinv_numerical", py::overload_cast<const cn&, const int>(&Lielab::functions::dexpinv_numerical<cn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", py::overload_cast<const glc&, const int>(&Lielab::functions::dexpinv_numerical<glc>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", py::overload_cast<const glr&, const int>(&Lielab::functions::dexpinv_numerical<glr>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", py::overload_cast<const rn&, const int>(&Lielab::functions::dexpinv_numerical<rn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", py::overload_cast<const se&, const int>(&Lielab::functions::dexpinv_numerical<se>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", py::overload_cast<const so&, const int>(&Lielab::functions::dexpinv_numerical<so>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", py::overload_cast<const sp&, const int>(&Lielab::functions::dexpinv_numerical<sp>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", py::overload_cast<const su&, const int>(&Lielab::functions::dexpinv_numerical<su>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::dexpinv_numerical<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("order") = 5);
    
    m_functions.def("dexpinv", py::overload_cast<const cn&, const int>(&Lielab::functions::dexpinv<cn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", py::overload_cast<const glc&, const int>(&Lielab::functions::dexpinv<glc>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", py::overload_cast<const glr&, const int>(&Lielab::functions::dexpinv<glr>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", py::overload_cast<const rn&, const int>(&Lielab::functions::dexpinv<rn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", py::overload_cast<const se&, const int>(&Lielab::functions::dexpinv<se>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", py::overload_cast<const so&, const int>(&Lielab::functions::dexpinv<so>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", py::overload_cast<const sp&, const int>(&Lielab::functions::dexpinv<sp>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", py::overload_cast<const su&, const int>(&Lielab::functions::dexpinv<su>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::dexpinv<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("order") = 5);

    m_functions.def("dlog_numerical", py::overload_cast<const cn&, const int>(&Lielab::functions::dlog_numerical<cn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", py::overload_cast<const glc&, const int>(&Lielab::functions::dlog_numerical<glc>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", py::overload_cast<const glr&, const int>(&Lielab::functions::dlog_numerical<glr>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", py::overload_cast<const rn&, const int>(&Lielab::functions::dlog_numerical<rn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", py::overload_cast<const se&, const int>(&Lielab::functions::dlog_numerical<se>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", py::overload_cast<const so&, const int>(&Lielab::functions::dlog_numerical<so>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", py::overload_cast<const sp&, const int>(&Lielab::functions::dlog_numerical<sp>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", py::overload_cast<const su&, const int>(&Lielab::functions::dlog_numerical<su>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::dlog_numerical<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("order") = 5);
    
    m_functions.def("dlog", py::overload_cast<const cn&, const int>(&Lielab::functions::dlog<cn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", py::overload_cast<const glc&, const int>(&Lielab::functions::dlog<glc>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", py::overload_cast<const glr&, const int>(&Lielab::functions::dlog<glr>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", py::overload_cast<const rn&, const int>(&Lielab::functions::dlog<rn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", py::overload_cast<const se&, const int>(&Lielab::functions::dlog<se>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", py::overload_cast<const so&, const int>(&Lielab::functions::dlog<so>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", py::overload_cast<const sp&, const int>(&Lielab::functions::dlog<sp>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", py::overload_cast<const su&, const int>(&Lielab::functions::dlog<su>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::dlog<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("order") = 5);

    m_functions.def("dloginv_numerical", py::overload_cast<const cn&, const int>(&Lielab::functions::dloginv_numerical<cn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv_numerical", py::overload_cast<const glc&, const int>(&Lielab::functions::dloginv_numerical<glc>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv_numerical", py::overload_cast<const glr&, const int>(&Lielab::functions::dloginv_numerical<glr>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv_numerical", py::overload_cast<const rn&, const int>(&Lielab::functions::dloginv_numerical<rn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv_numerical", py::overload_cast<const se&, const int>(&Lielab::functions::dloginv_numerical<se>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv_numerical", py::overload_cast<const so&, const int>(&Lielab::functions::dloginv_numerical<so>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv_numerical", py::overload_cast<const sp&, const int>(&Lielab::functions::dloginv_numerical<sp>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv_numerical", py::overload_cast<const su&, const int>(&Lielab::functions::dloginv_numerical<su>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv_numerical", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::dloginv_numerical<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("order") = 5);
    
    m_functions.def("dloginv", py::overload_cast<const cn&, const int>(&Lielab::functions::dloginv<cn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv", py::overload_cast<const glc&, const int>(&Lielab::functions::dloginv<glc>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv", py::overload_cast<const glr&, const int>(&Lielab::functions::dloginv<glr>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv", py::overload_cast<const rn&, const int>(&Lielab::functions::dloginv<rn>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv", py::overload_cast<const se&, const int>(&Lielab::functions::dloginv<se>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv", py::overload_cast<const so&, const int>(&Lielab::functions::dloginv<so>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv", py::overload_cast<const sp&, const int>(&Lielab::functions::dloginv<sp>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv", py::overload_cast<const su&, const int>(&Lielab::functions::dloginv<su>), py::arg("a"), py::arg("order") = 5);
    m_functions.def("dloginv", py::overload_cast<const CompositeAlgebra&, const int>(&Lielab::functions::dloginv<CompositeAlgebra, CompositeAlgebra>), py::arg("a"), py::arg("order") = 5);
}
