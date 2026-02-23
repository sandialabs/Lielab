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

#include "pybind_GLC.hpp"

namespace py = pybind11;

void bind_GLC(py::module& m_domain)
{
    auto Lielab_domain_GLC = py::class_<Lielab::domain::GLC>(m_domain, "GLC");
    Lielab_domain_GLC.def_readwrite("point", &Lielab::domain::GLC::point);
    Lielab_domain_GLC.def(py::init<>());
    Lielab_domain_GLC.def(py::init<const Eigen::MatrixXcd&>());
    Lielab_domain_GLC.def_static("identity", &Lielab::domain::GLC::identity);
    Lielab_domain_GLC.def_static("project", &Lielab::domain::GLC::project);
    Lielab_domain_GLC.def(py::init<const int>());
    Lielab_domain_GLC.def("to_string", &Lielab::domain::GLC::to_string);
    Lielab_domain_GLC.def("__repr__", [](const Lielab::domain::GLC& self)
        {
            return "<lielab.domain.GLC>";
        });
    Lielab_domain_GLC.def("__str__", [](const Lielab::domain::GLC& self)
        {
            return "<lielab.domain.GLC>";
        });
    Lielab_domain_GLC.def("get_dimension", &Lielab::domain::GLC::get_dimension);
    Lielab_domain_GLC.def("get_size", &Lielab::domain::GLC::get_size);
    Lielab_domain_GLC.def("is_abelian", &Lielab::domain::GLC::is_abelian);
    Lielab_domain_GLC.def("get_shape", &Lielab::domain::GLC::get_shape);
    Lielab_domain_GLC.def("get_point", &Lielab::domain::GLC::get_point);
    Lielab_domain_GLC.def("serialize", &Lielab::domain::GLC::serialize);
    Lielab_domain_GLC.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::GLC::unserialize));
    Lielab_domain_GLC.def("get_matrix", &Lielab::domain::GLC::get_matrix);
    Lielab_domain_GLC.def("__call__", [](const Lielab::domain::GLC& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_GLC.def(py::self * py::self);
    Lielab_domain_GLC.def(py::self *= py::self);
    Lielab_domain_GLC.def("inverse", &Lielab::domain::GLC::inverse);

    // Other python methods
    Lielab_domain_GLC.def(py::pickle([](const Lielab::domain::GLC& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 2)
                throw std::runtime_error("GLC: Invalid state.");

            Lielab::domain::GLC obj;
            obj.point = t[0].cast<Lielab::domain::GLC::point_t>();
            return obj;
        }));
}
