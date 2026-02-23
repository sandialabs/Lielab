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

#include "pybind_SE.hpp"

namespace py = pybind11;

void bind_SE(py::module& m_domain)
{
    auto Lielab_domain_SE = py::class_<Lielab::domain::SE>(m_domain, "SE");
    Lielab_domain_SE.def_readwrite("point", &Lielab::domain::SE::point);
    Lielab_domain_SE.def(py::init<>());
    Lielab_domain_SE.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_SE.def_static("identity", &Lielab::domain::SE::identity);
    Lielab_domain_SE.def_static("project", &Lielab::domain::SE::project);
    Lielab_domain_SE.def(py::init<const int>());
    Lielab_domain_SE.def(py::init<const Lielab::domain::RN&, const Lielab::domain::SO&>());
    Lielab_domain_SE.def("to_string", &Lielab::domain::SE::to_string);
    Lielab_domain_SE.def("__repr__", [](const Lielab::domain::SE& self)
        {
            return "<lielab.domain.SE>";
        });
    Lielab_domain_SE.def("__str__", [](const Lielab::domain::SE& self)
        {
            return "<lielab.domain.SE>";
        });
    Lielab_domain_SE.def("get_dimension", &Lielab::domain::SE::get_dimension);
    Lielab_domain_SE.def("get_size", &Lielab::domain::SE::get_size);
    Lielab_domain_SE.def("is_abelian", &Lielab::domain::SE::is_abelian);
    Lielab_domain_SE.def("get_shape", &Lielab::domain::SE::get_shape);
    Lielab_domain_SE.def("get_point", &Lielab::domain::SE::get_point);
    Lielab_domain_SE.def("get_data", &Lielab::domain::SE::get_point);
    Lielab_domain_SE.def("serialize", &Lielab::domain::SE::serialize);
    Lielab_domain_SE.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::SE::unserialize));
    Lielab_domain_SE.def("get_matrix", &Lielab::domain::SE::get_matrix);
    Lielab_domain_SE.def("__call__", [](const Lielab::domain::SE& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SE.def(py::self * py::self);
    Lielab_domain_SE.def(py::self *= py::self);
    Lielab_domain_SE.def("inverse", &Lielab::domain::SE::inverse);

    // Misc python functions
    Lielab_domain_SE.def(py::pickle([](const Lielab::domain::SE& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point, obj._shape);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 2)
                throw std::runtime_error("SE: Invalid state.");

            Lielab::domain::SE obj;
            obj.point = t[0].cast<Lielab::domain::SE::point_t>();
            obj._shape = t[1].cast<int>();
            return obj;
        }));
}
