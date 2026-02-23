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

#include "pybind_so.hpp"

namespace py = pybind11;

void bind_so(py::module& m_domain)
{
    /*!
    * Bindings for Lielab::domain::so
    */

    auto Lielab_domain_so = py::class_<Lielab::domain::so>(m_domain, "so");
    Lielab_domain_so.def_readwrite("point", &Lielab::domain::so::point);
    Lielab_domain_so.def(py::init<>());
    Lielab_domain_so.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_so.def("basis", &Lielab::domain::so::basis);
    Lielab_domain_so.def("zero", &Lielab::domain::so::zero);
    Lielab_domain_so.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::so::from_vector));
    Lielab_domain_so.def("project", &Lielab::domain::so::project);
    Lielab_domain_so.def(py::init<const int>());
    Lielab_domain_so.def("to_string", &Lielab::domain::so::to_string);
    Lielab_domain_so.def("__repr__", [](const Lielab::domain::so& self)
        {
            return "<lielab.domain.so>";
        });
    Lielab_domain_so.def("__str__", [](const Lielab::domain::so& self)
        {
            return "<lielab.domain.so>";
        });
    Lielab_domain_so.def("get_dimension", &Lielab::domain::so::get_dimension);
    Lielab_domain_so.def("get_size", &Lielab::domain::so::get_size);
    Lielab_domain_so.def("is_abelian", &Lielab::domain::so::is_abelian);
    Lielab_domain_so.def("get_shape", &Lielab::domain::so::get_shape);
    Lielab_domain_so.def("get_point", &Lielab::domain::so::get_point);
    Lielab_domain_so.def("serialize", &Lielab::domain::so::serialize);
    Lielab_domain_so.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::so::unserialize));
    Lielab_domain_so.def("get_matrix", &Lielab::domain::so::get_matrix);
    Lielab_domain_so.def("get_vector", &Lielab::domain::so::get_vector);
    Lielab_domain_so.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::so::set_vector));
    Lielab_domain_so.def("__call__", [](const Lielab::domain::so& self, const int index)
        {
            return self(index);
        });
    Lielab_domain_so.def("__call__", [](const Lielab::domain::so& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_so.def(py::self + py::self);
    Lielab_domain_so.def(py::self += py::self);
    Lielab_domain_so.def(py::self - py::self);
    Lielab_domain_so.def(py::self -= py::self);
    Lielab_domain_so.def(-py::self);
    Lielab_domain_so.def(py::self * int());
    Lielab_domain_so.def(py::self * double());
    Lielab_domain_so.def(int() * py::self);
    Lielab_domain_so.def(double() * py::self);
    Lielab_domain_so.def(py::self *= int());
    Lielab_domain_so.def(py::self *= double());
    Lielab_domain_so.def(py::self / int());
    Lielab_domain_so.def(py::self / double());
    Lielab_domain_so.def(py::self /= int());
    Lielab_domain_so.def(py::self /= double());

    // Other misc python methods
    Lielab_domain_so.def(py::pickle([](const Lielab::domain::so& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 1)
                throw std::runtime_error("so: Invalid state.");

            Lielab::domain::so obj;
            obj.point = t[0].cast<Lielab::domain::so::point_t>();
            return obj;
        }));
}
