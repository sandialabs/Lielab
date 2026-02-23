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

#include "pybind_SU.hpp"

namespace py = pybind11;

void bind_SU(py::module& m_domain)
{
    auto Lielab_domain_SU = py::class_<Lielab::domain::SU>(m_domain, "SU");
    Lielab_domain_SU.def_readwrite("point", &Lielab::domain::SU::point);
    Lielab_domain_SU.def(py::init<>());
    Lielab_domain_SU.def(py::init<const Eigen::MatrixXcd&>());
    Lielab_domain_SU.def_static("identity", Lielab::domain::SU::identity);
    // .def_static("project", &Lielab::domain::SU::project); // TODO:
    Lielab_domain_SU.def(py::init<const int>());
    Lielab_domain_SU.def_static("from_quaternion", &Lielab::domain::SU::from_quaternion);
    Lielab_domain_SU.def_static("from_SO3", &Lielab::domain::SU::from_SO3);
    Lielab_domain_SU.def("to_string", &Lielab::domain::SU::to_string);
    Lielab_domain_SU.def("__repr__", [](const Lielab::domain::SU & self)
        {
            return "<lielab.domain.SU>";
        });
    Lielab_domain_SU.def("__str__", [](const Lielab::domain::SU & self)
        {
            return "<lielab.domain.SU>";
        });
    Lielab_domain_SU.def("get_dimension", &Lielab::domain::SU::get_dimension);
    Lielab_domain_SU.def("get_size", &Lielab::domain::SU::get_size);
    Lielab_domain_SU.def("is_abelian", &Lielab::domain::SU::is_abelian);
    Lielab_domain_SU.def("get_shape", &Lielab::domain::SU::get_shape);
    Lielab_domain_SU.def("get_point", &Lielab::domain::SU::get_point);
    Lielab_domain_SU.def("serialize", &Lielab::domain::SU::serialize);
    Lielab_domain_SU.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::SU::unserialize));
    Lielab_domain_SU.def("get_matrix", &Lielab::domain::SU::get_matrix);
    Lielab_domain_SU.def("__call__", [](const Lielab::domain::SU& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SU.def("to_quaternion", &Lielab::domain::SU::to_quaternion);
    Lielab_domain_SU.def(py::self * py::self);
    Lielab_domain_SU.def(py::self *= py::self);
    Lielab_domain_SU.def("inverse", &Lielab::domain::SU::inverse);

    // Other misc python
    Lielab_domain_SU.def(py::pickle([](const Lielab::domain::SU& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 1)
                throw std::runtime_error("SU: Invalid state.");

            Lielab::domain::SU obj;
            obj.point = t[0].cast<Lielab::domain::SU::point_t>();
            return obj;
        }));
}
