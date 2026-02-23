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

#include "pybind_glr.hpp"

namespace py = pybind11;

void bind_glr(py::module& m_domain)
{
    /*!
    * Bindings for Lielab::domain::glr
    */

    auto Lielab_domain_glr = py::class_<Lielab::domain::glr>(m_domain, "glr");
    Lielab_domain_glr.def_readwrite("point", &Lielab::domain::glr::point);
    Lielab_domain_glr.def(py::init<>());
    Lielab_domain_glr.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_glr.def_static("basis", &Lielab::domain::glr::basis);
    Lielab_domain_glr.def_static("zero", &Lielab::domain::glr::zero);
    Lielab_domain_glr.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::glr::from_vector));
    Lielab_domain_glr.def_static("project", &Lielab::domain::glr::project);
    Lielab_domain_glr.def(py::init<const int>());
    Lielab_domain_glr.def("to_string", &Lielab::domain::glr::to_string);
    Lielab_domain_glr.def("__repr__", [](const Lielab::domain::glr& self)
        {
            return "<lielab.domain.glr>";
        });
    Lielab_domain_glr.def("__str__", [](const Lielab::domain::glr& self)
        {
            return "<lielab.domain.glr>";
        });
    Lielab_domain_glr.def("get_dimension", &Lielab::domain::glr::get_dimension);
    Lielab_domain_glr.def("get_size", &Lielab::domain::glr::get_size);
    Lielab_domain_glr.def("is_abelian", &Lielab::domain::glr::is_abelian);
    Lielab_domain_glr.def("get_shape", &Lielab::domain::glr::get_shape);
    Lielab_domain_glr.def("get_point", &Lielab::domain::glr::get_point);
    Lielab_domain_glr.def("serialize", &Lielab::domain::glr::serialize);
    Lielab_domain_glr.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::glr::unserialize));
    Lielab_domain_glr.def("get_matrix", &Lielab::domain::glr::get_matrix);
    Lielab_domain_glr.def("get_vector", &Lielab::domain::glr::get_vector);
    Lielab_domain_glr.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::glr::set_vector));
    Lielab_domain_glr.def("__call__", [](const Lielab::domain::glr& self, const int index)
        {
            return self(index);
        });
    Lielab_domain_glr.def("__call__", [](const Lielab::domain::glr& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_glr.def(py::self + py::self);
    Lielab_domain_glr.def(py::self += py::self);
    Lielab_domain_glr.def(py::self - py::self);
    Lielab_domain_glr.def(py::self -= py::self);
    Lielab_domain_glr.def(-py::self);
    Lielab_domain_glr.def(py::self * int());
    Lielab_domain_glr.def(py::self * double());
    Lielab_domain_glr.def(py::self *= int());
    Lielab_domain_glr.def(py::self *= double());
    Lielab_domain_glr.def(py::self / int());
    Lielab_domain_glr.def(py::self / double());
    Lielab_domain_glr.def(py::self /= int());
    Lielab_domain_glr.def(py::self /= double());
    Lielab_domain_glr.def(int() * py::self);
    Lielab_domain_glr.def(double() * py::self);

    // Other misc python
    Lielab_domain_glr.def(py::pickle([](const Lielab::domain::glr& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 1)
                throw std::runtime_error("glr: Invalid state.");

            Lielab::domain::glr obj;
            obj.point = t[0].cast<Lielab::domain::glr::point_t>();
            return obj;
        }));
}
