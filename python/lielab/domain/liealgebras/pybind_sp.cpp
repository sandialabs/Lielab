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

#include "pybind_sp.hpp"

namespace py = pybind11;

void bind_sp(py::module& m_domain)
{
    /*!
    * Bindings for Lielab::domain::sp
    */
    
    auto Lielab_domain_sp = py::class_<Lielab::domain::sp>(m_domain, "sp");
    Lielab_domain_sp.def_readwrite("point", &Lielab::domain::sp::point);
    Lielab_domain_sp.def(py::init<>());
    Lielab_domain_sp.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_sp.def_static("basis", &Lielab::domain::sp::basis);
    Lielab_domain_sp.def_static("zero", &Lielab::domain::sp::zero);
    Lielab_domain_sp.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::sp::from_vector));
    Lielab_domain_sp.def_static("project", &Lielab::domain::sp::project);
    Lielab_domain_sp.def(py::init<const int>());
    Lielab_domain_sp.def("to_string", &Lielab::domain::sp::to_string);
    Lielab_domain_sp.def("__repr__", [](const Lielab::domain::sp& self)
        {
            return "<lielab.domain.sp>";
        });
    Lielab_domain_sp.def("__str__", [](const Lielab::domain::sp& self)
        {
            return "<lielab.domain.sp>";
        });
    Lielab_domain_sp.def("get_dimension", &Lielab::domain::sp::get_dimension);
    Lielab_domain_sp.def("get_size", &Lielab::domain::sp::get_size);
    Lielab_domain_sp.def("is_abelian", &Lielab::domain::sp::is_abelian);
    Lielab_domain_sp.def("get_shape", &Lielab::domain::sp::get_shape);
    Lielab_domain_sp.def("get_point", &Lielab::domain::sp::get_point);
    Lielab_domain_sp.def("serialize", &Lielab::domain::sp::serialize);
    Lielab_domain_sp.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::sp::unserialize));
    Lielab_domain_sp.def("get_matrix", &Lielab::domain::sp::get_matrix);
    Lielab_domain_sp.def("get_vector", &Lielab::domain::sp::get_vector);
    Lielab_domain_sp.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::sp::set_vector));
    Lielab_domain_sp.def("__call__", [](const Lielab::domain::sp& self, const int index)
        {
            return self(index);
        });
    Lielab_domain_sp.def("__call__", [](const Lielab::domain::sp& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_sp.def(py::self + py::self);
    Lielab_domain_sp.def(py::self += py::self);
    Lielab_domain_sp.def(py::self - py::self);
    Lielab_domain_sp.def(py::self -= py::self);
    Lielab_domain_sp.def(-py::self);
    Lielab_domain_sp.def(py::self * int());
    Lielab_domain_sp.def(py::self * double());
    Lielab_domain_sp.def(int() * py::self);
    Lielab_domain_sp.def(double() * py::self);
    Lielab_domain_sp.def(py::self *= int());
    Lielab_domain_sp.def(py::self *= double());
    Lielab_domain_sp.def(py::self / int());
    Lielab_domain_sp.def(py::self / double());
    Lielab_domain_sp.def(py::self /= int());
    Lielab_domain_sp.def(py::self /= double());

    // Other misc python methods
    Lielab_domain_sp.def(py::pickle([](const Lielab::domain::sp& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 1)
                throw std::runtime_error("sp: Invalid state.");

            Lielab::domain::sp obj;
            obj.point = t[0].cast<Lielab::domain::sp::point_t>();
            return obj;
        }));
}
