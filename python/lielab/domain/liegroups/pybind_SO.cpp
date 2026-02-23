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

#include "pybind_SO.hpp"

namespace py = pybind11;

void bind_SO(py::module& m_domain)
{
    auto Lielab_domain_SO = py::class_<Lielab::domain::SO>(m_domain, "SO");
    Lielab_domain_SO.def_readwrite("point", &Lielab::domain::SO::point);
    Lielab_domain_SO.def(py::init<>());
    Lielab_domain_SO.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_SO.def_static("identity", &Lielab::domain::SO::identity);
    Lielab_domain_SO.def_static("project", &Lielab::domain::SO::project);
    Lielab_domain_SO.def(py::init<const int>());
    Lielab_domain_SO.def_static("from_eulerangles_body123", &Lielab::domain::SO::from_eulerangles_body123, "The from_eulerangles_body123 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body231", &Lielab::domain::SO::from_eulerangles_body231, "The from_eulerangles_body231 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body312", &Lielab::domain::SO::from_eulerangles_body312, "The from_eulerangles_body312 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body132", &Lielab::domain::SO::from_eulerangles_body132, "The from_eulerangles_body132 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body213", &Lielab::domain::SO::from_eulerangles_body213, "The from_eulerangles_body213 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body321", &Lielab::domain::SO::from_eulerangles_body321, "The from_eulerangles_body321 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body121", &Lielab::domain::SO::from_eulerangles_body121, "The from_eulerangles_body121 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body131", &Lielab::domain::SO::from_eulerangles_body131, "The from_eulerangles_body131 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body212", &Lielab::domain::SO::from_eulerangles_body212, "The from_eulerangles_body212 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body232", &Lielab::domain::SO::from_eulerangles_body232, "The from_eulerangles_body232 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body313", &Lielab::domain::SO::from_eulerangles_body313, "The from_eulerangles_body313 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body323", &Lielab::domain::SO::from_eulerangles_body323, "The from_eulerangles_body323 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space123", &Lielab::domain::SO::from_eulerangles_space123, "The from_eulerangles_space123 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space231", &Lielab::domain::SO::from_eulerangles_space231, "The from_eulerangles_space231 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space312", &Lielab::domain::SO::from_eulerangles_space312, "The from_eulerangles_space312 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space132", &Lielab::domain::SO::from_eulerangles_space132, "The from_eulerangles_space132 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space213", &Lielab::domain::SO::from_eulerangles_space213, "The from_eulerangles_space213 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space321", &Lielab::domain::SO::from_eulerangles_space321, "The from_eulerangles_space321 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space121", &Lielab::domain::SO::from_eulerangles_space121, "The from_eulerangles_space121 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space131", &Lielab::domain::SO::from_eulerangles_space131, "The from_eulerangles_space131 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space212", &Lielab::domain::SO::from_eulerangles_space212, "The from_eulerangles_space212 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space232", &Lielab::domain::SO::from_eulerangles_space232, "The from_eulerangles_space232 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space313", &Lielab::domain::SO::from_eulerangles_space313, "The from_eulerangles_space313 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space323", &Lielab::domain::SO::from_eulerangles_space323, "The from_eulerangles_space323 function.");
    Lielab_domain_SO.def_static("from_quaternion", &Lielab::domain::SO::from_quaternion, "The from_quaternion function.");
    Lielab_domain_SO.def_static("from_rodriguesvector", &Lielab::domain::SO::from_rodriguesvector, "The from_rodriguesvector function.");
    Lielab_domain_SO.def_static("from_SU2", &Lielab::domain::SO::from_SU2);
    Lielab_domain_SO.def("to_string", &Lielab::domain::SO::to_string);
    Lielab_domain_SO.def("__repr__", [](const Lielab::domain::SO& self)
        {
            return "<lielab.domain.SO>";
        });
    Lielab_domain_SO.def("__str__", [](const Lielab::domain::SO& self)
        {
            return "<lielab.domain.SO>";
        });
    Lielab_domain_SO.def("get_dimension", &Lielab::domain::SO::get_dimension);
    Lielab_domain_SO.def("get_size", &Lielab::domain::SO::get_size);
    Lielab_domain_SO.def("is_abelian", &Lielab::domain::SO::is_abelian);
    Lielab_domain_SO.def("get_shape", &Lielab::domain::SO::get_shape);
    Lielab_domain_SO.def("get_point", &Lielab::domain::SO::get_point);
    Lielab_domain_SO.def("serialize", &Lielab::domain::SO::serialize);
    Lielab_domain_SO.def("unserialize", py::overload_cast<const Eigen::VectorXd &>(&Lielab::domain::SO::unserialize));
    Lielab_domain_SO.def("get_matrix", &Lielab::domain::SO::get_matrix);
    Lielab_domain_SO.def("__call__", [](const Lielab::domain::SO & self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SO.def("to_eulerangles_body123", &Lielab::domain::SO::to_eulerangles_body123, "The to_eulerangles_body123 function.");
    Lielab_domain_SO.def("to_eulerangles_body231", &Lielab::domain::SO::to_eulerangles_body231, "The to_eulerangles_body231 function.");
    Lielab_domain_SO.def("to_eulerangles_body312", &Lielab::domain::SO::to_eulerangles_body312, "The to_eulerangles_body312 function.");
    Lielab_domain_SO.def("to_eulerangles_body132", &Lielab::domain::SO::to_eulerangles_body132, "The to_eulerangles_body132 function.");
    Lielab_domain_SO.def("to_eulerangles_body213", &Lielab::domain::SO::to_eulerangles_body213, "The to_eulerangles_body213 function.");
    Lielab_domain_SO.def("to_eulerangles_body321", &Lielab::domain::SO::to_eulerangles_body321, "The to_eulerangles_body321 function.");
    Lielab_domain_SO.def("to_eulerangles_body121", &Lielab::domain::SO::to_eulerangles_body121, "The to_eulerangles_body121 function.");
    Lielab_domain_SO.def("to_eulerangles_body131", &Lielab::domain::SO::to_eulerangles_body131, "The to_eulerangles_body131 function.");
    Lielab_domain_SO.def("to_eulerangles_body212", &Lielab::domain::SO::to_eulerangles_body212, "The to_eulerangles_body212 function.");
    Lielab_domain_SO.def("to_eulerangles_body232", &Lielab::domain::SO::to_eulerangles_body232, "The to_eulerangles_body232 function.");
    Lielab_domain_SO.def("to_eulerangles_body313", &Lielab::domain::SO::to_eulerangles_body313, "The to_eulerangles_body313 function.");
    Lielab_domain_SO.def("to_eulerangles_body323", &Lielab::domain::SO::to_eulerangles_body323, "The to_eulerangles_body323 function.");
    Lielab_domain_SO.def("to_eulerangles_space123", &Lielab::domain::SO::to_eulerangles_space123, "The to_eulerangles_space123 function.");
    Lielab_domain_SO.def("to_eulerangles_space231", &Lielab::domain::SO::to_eulerangles_space231, "The to_eulerangles_space231 function.");
    Lielab_domain_SO.def("to_eulerangles_space312", &Lielab::domain::SO::to_eulerangles_space312, "The to_eulerangles_space312 function.");
    Lielab_domain_SO.def("to_eulerangles_space132", &Lielab::domain::SO::to_eulerangles_space132, "The to_eulerangles_space132 function.");
    Lielab_domain_SO.def("to_eulerangles_space213", &Lielab::domain::SO::to_eulerangles_space213, "The to_eulerangles_space213 function.");
    Lielab_domain_SO.def("to_eulerangles_space321", &Lielab::domain::SO::to_eulerangles_space321, "The to_eulerangles_space321 function.");
    Lielab_domain_SO.def("to_eulerangles_space121", &Lielab::domain::SO::to_eulerangles_space121, "The to_eulerangles_space121 function.");
    Lielab_domain_SO.def("to_eulerangles_space131", &Lielab::domain::SO::to_eulerangles_space131, "The to_eulerangles_space131 function.");
    Lielab_domain_SO.def("to_eulerangles_space212", &Lielab::domain::SO::to_eulerangles_space212, "The to_eulerangles_space212 function.");
    Lielab_domain_SO.def("to_eulerangles_space232", &Lielab::domain::SO::to_eulerangles_space232, "The to_eulerangles_space232 function.");
    Lielab_domain_SO.def("to_eulerangles_space313", &Lielab::domain::SO::to_eulerangles_space313, "The to_eulerangles_space313 function.");
    Lielab_domain_SO.def("to_eulerangles_space323", &Lielab::domain::SO::to_eulerangles_space323, "The to_eulerangles_space323 function.");
    Lielab_domain_SO.def("to_quaternion", &Lielab::domain::SO::to_quaternion);
    Lielab_domain_SO.def("to_gibbs", &Lielab::domain::SO::to_gibbs);
    Lielab_domain_SO.def(py::self * py::self);
    Lielab_domain_SO.def(py::self *= py::self);
    Lielab_domain_SO.def("inverse", &Lielab::domain::SO::inverse);
    Lielab_domain_SO.def("project_onto_tangent_space", &Lielab::domain::SO::project_onto_tangent_space);

    // Other misc functions
    Lielab_domain_SO.def(py::pickle([](const Lielab::domain::SO& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 1)
                throw std::runtime_error("SO: Invalid state.");

            Lielab::domain::SO obj;
            obj.point = t[0].cast<Lielab::domain::SO::point_t>();
            return obj;
        }));
}
