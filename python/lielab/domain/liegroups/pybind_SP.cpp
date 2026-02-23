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

#include "pybind_SP.hpp"

namespace py = pybind11;

void bind_SP(py::module& m_domain)
{
    auto Lielab_domain_SP = py::class_<Lielab::domain::SP>(m_domain, "SP");
    Lielab_domain_SP.def_readwrite("point", &Lielab::domain::SP::point);
    Lielab_domain_SP.def(py::init<>());
    Lielab_domain_SP.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_SP.def_static("identity", Lielab::domain::SP::identity);
    // .def_static("project", &Lielab::domain::SP::project); // TODO:
    Lielab_domain_SP.def(py::init<const int>());
    Lielab_domain_SP.def("to_string", &Lielab::domain::SP::to_string);
    Lielab_domain_SP.def("__repr__", [](const Lielab::domain::SP& self)
        {
            return "<lielab.domain.SP>";
        });
    Lielab_domain_SP.def("__str__", [](const Lielab::domain::SP& self)
        {
            return "<lielab.domain.SP>";
        });
    Lielab_domain_SP.def("get_dimension", &Lielab::domain::SP::get_dimension);
    Lielab_domain_SP.def("get_size", &Lielab::domain::SP::get_size);
    Lielab_domain_SP.def("is_abelian", &Lielab::domain::SP::is_abelian);
    Lielab_domain_SP.def("get_shape", &Lielab::domain::SP::get_shape);
    Lielab_domain_SP.def("get_point", &Lielab::domain::SP::get_point);
    Lielab_domain_SP.def("serialize", &Lielab::domain::SP::serialize);
    Lielab_domain_SP.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::SP::unserialize));
    Lielab_domain_SP.def("get_matrix", &Lielab::domain::SP::get_matrix);
    Lielab_domain_SP.def("__call__", [](const Lielab::domain::SP& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SP.def(py::self * py::self);
    Lielab_domain_SP.def(py::self *= py::self);
    Lielab_domain_SP.def("inverse", &Lielab::domain::SP::inverse);

    // Other python funcs
    Lielab_domain_SP.def(py::pickle([](const Lielab::domain::SP& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 1)
                throw std::runtime_error("SP: Invalid state.");

            Lielab::domain::SP obj;
            obj.point = t[0].cast<Lielab::domain::SP::point_t>();
            return obj;
        }));
}
