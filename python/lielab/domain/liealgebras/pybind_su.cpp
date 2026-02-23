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

#include "pybind_su.hpp"

namespace py = pybind11;

void bind_su(py::module& m_domain)
{
    /*!
    * Bindings for Lielab::domain::su
    */
    
    auto Lielab_domain_su = py::class_<Lielab::domain::su>(m_domain, "su");
    Lielab_domain_su.def_readwrite("point", &Lielab::domain::su::point);
    Lielab_domain_su.def(py::init<>());
    Lielab_domain_su.def(py::init<const Eigen::MatrixXcd&>());
    Lielab_domain_su.def_static("basis", &Lielab::domain::su::basis);
    Lielab_domain_su.def_static("zero", &Lielab::domain::su::zero);
    Lielab_domain_su.def_static("project", &Lielab::domain::su::project);
    Lielab_domain_su.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::su::from_vector));
    Lielab_domain_su.def(py::init<const int>());
    Lielab_domain_su.def("to_string", &Lielab::domain::su::to_string);
    Lielab_domain_su.def("__repr__", [](const Lielab::domain::su& self)
        {
            return "<lielab.domain.su>";
        });
    Lielab_domain_su.def("__str__", [](const Lielab::domain::su& self)
        {
            return "<lielab.domain.su>";
        });
    Lielab_domain_su.def("get_dimension", &Lielab::domain::su::get_dimension);
    Lielab_domain_su.def("get_size", &Lielab::domain::su::get_size);
    Lielab_domain_su.def("is_abelian", &Lielab::domain::su::is_abelian);
    Lielab_domain_su.def("get_shape", &Lielab::domain::su::get_shape);
    Lielab_domain_su.def("get_point", &Lielab::domain::su::get_point);
    Lielab_domain_su.def("serialize", &Lielab::domain::su::serialize);
    Lielab_domain_su.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::su::unserialize));
    Lielab_domain_su.def("get_matrix", &Lielab::domain::su::get_matrix);
    Lielab_domain_su.def("get_vector", &Lielab::domain::su::get_vector);
    Lielab_domain_su.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::su::set_vector));
    Lielab_domain_su.def("__call__", [](const Lielab::domain::su& self, const int index)
        {
            return self(index);
        });
    Lielab_domain_su.def("__call__", [](const Lielab::domain::su& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_su.def(py::self + py::self);
    Lielab_domain_su.def(py::self += py::self);
    Lielab_domain_su.def(py::self - py::self);
    Lielab_domain_su.def(py::self -= py::self);
    Lielab_domain_su.def(-py::self);
    Lielab_domain_su.def(py::self * int());
    Lielab_domain_su.def(py::self * double());
    Lielab_domain_su.def(int() * py::self);
    Lielab_domain_su.def(double() * py::self);
    Lielab_domain_su.def(py::self *= int());
    Lielab_domain_su.def(py::self *= double());
    Lielab_domain_su.def(py::self / int());
    Lielab_domain_su.def(py::self / double());
    Lielab_domain_su.def(py::self /= int());
    Lielab_domain_su.def(py::self /= double());

    // other python methods
    Lielab_domain_su.def(py::pickle([](const Lielab::domain::su& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 1)
                throw std::runtime_error("su: Invalid state.");

            Lielab::domain::su obj;
            obj.point = t[0].cast<Lielab::domain::su::point_t>();
            return obj;
        }));
}
