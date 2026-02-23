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

#include "pybind_GLR.hpp"

namespace py = pybind11;

void bind_GLR(py::module& m_domain)
{
    auto Lielab_domain_GLR = py::class_<Lielab::domain::GLR>(m_domain, "GLR");
    Lielab_domain_GLR.def_readwrite("point", &Lielab::domain::GLR::point);
    Lielab_domain_GLR.def(py::init<>());
    Lielab_domain_GLR.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_GLR.def_static("identity", &Lielab::domain::GLR::identity);
    Lielab_domain_GLR.def_static("project", &Lielab::domain::GLR::project);
    Lielab_domain_GLR.def(py::init<const int>());
    Lielab_domain_GLR.def("to_string", &Lielab::domain::GLR::to_string);
    Lielab_domain_GLR.def("__repr__", [](const Lielab::domain::GLR& self)
        {
            return "<lielab.domain.GLR>";
        });
    Lielab_domain_GLR.def("__str__", [](const Lielab::domain::GLR& self)
        {
            return "<lielab.domain.GLR>";
        });
    Lielab_domain_GLR.def("get_dimension", &Lielab::domain::GLR::get_dimension);
    Lielab_domain_GLR.def("get_size", &Lielab::domain::GLR::get_size);
    Lielab_domain_GLR.def("is_abelian", &Lielab::domain::GLR::is_abelian);
    Lielab_domain_GLR.def("get_shape", &Lielab::domain::GLR::get_shape);
    Lielab_domain_GLR.def("get_point", &Lielab::domain::GLR::get_point);
    Lielab_domain_GLR.def("serialize", &Lielab::domain::GLR::serialize);
    Lielab_domain_GLR.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::GLR::unserialize));
    Lielab_domain_GLR.def("get_matrix", &Lielab::domain::GLR::get_matrix);
    Lielab_domain_GLR.def("__call__", [](const Lielab::domain::GLR& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_GLR.def(py::self * py::self);
    Lielab_domain_GLR.def(py::self *= py::self);
    Lielab_domain_GLR.def("inverse", &Lielab::domain::GLR::inverse);
    
    // Misc python methods
    Lielab_domain_GLR.def(py::pickle([](const Lielab::domain::GLR& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 1)
                throw std::runtime_error("GLR: Invalid state.");

            Lielab::domain::GLR obj;
            obj.point = t[0].cast<Lielab::domain::GLR::point_t>();
            return obj;
        }));
}
