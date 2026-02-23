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

#include "pybind_se.hpp"

namespace py = pybind11;

void bind_se(py::module& m_domain)
{
    /*!
    * Bindings for Lielab::domain::se
    */

    auto Lielab_domain_se = py::class_<Lielab::domain::se>(m_domain, "se");
    Lielab_domain_se.def_readwrite("point", &Lielab::domain::se::point);
    Lielab_domain_se.def(py::init<>());
    Lielab_domain_se.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_se.def_static("basis", &Lielab::domain::se::basis);
    Lielab_domain_se.def_static("zero", &Lielab::domain::se::zero);
    Lielab_domain_se.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::se::from_vector));
    Lielab_domain_se.def_static("project", &Lielab::domain::se::project);
    Lielab_domain_se.def(py::init<const int>());
    Lielab_domain_se.def(py::init<const Lielab::domain::rn&, const Lielab::domain::so&>());
    Lielab_domain_se.def("to_string", &Lielab::domain::se::to_string);
    Lielab_domain_se.def("__repr__", [](const Lielab::domain::se& self)
        {
            return "<lielab.domain.se>";
        });
    Lielab_domain_se.def("__str__", [](const Lielab::domain::se& self)
        {
            return "<lielab.domain.se>";
        });
    Lielab_domain_se.def("get_dimension", &Lielab::domain::se::get_dimension);
    Lielab_domain_se.def("get_size", &Lielab::domain::se::get_size);
    Lielab_domain_se.def("is_abelian", &Lielab::domain::se::is_abelian);
    Lielab_domain_se.def("get_shape", &Lielab::domain::se::get_shape);
    Lielab_domain_se.def("get_point", &Lielab::domain::se::get_point);
    Lielab_domain_se.def("serialize", &Lielab::domain::se::serialize);
    Lielab_domain_se.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::se::unserialize));
    Lielab_domain_se.def("get_matrix", &Lielab::domain::se::get_matrix);
    Lielab_domain_se.def("get_vector", &Lielab::domain::se::get_vector);
    Lielab_domain_se.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::se::set_vector));
    Lielab_domain_se.def("__call__", [](const Lielab::domain::se& self, const int index)
        {
            return self(index);
        });
    Lielab_domain_se.def("__call__", [](const Lielab::domain::se& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_se.def(py::self + py::self);
    Lielab_domain_se.def(py::self += py::self);
    Lielab_domain_se.def(py::self - py::self);
    Lielab_domain_se.def(py::self -= py::self);
    Lielab_domain_se.def(-py::self);
    Lielab_domain_se.def(py::self * int());
    Lielab_domain_se.def(py::self * double());
    Lielab_domain_se.def(int() * py::self);
    Lielab_domain_se.def(double() * py::self);
    Lielab_domain_se.def(py::self *= int());
    Lielab_domain_se.def(py::self *= double());
    Lielab_domain_se.def(py::self / int());
    Lielab_domain_se.def(py::self / double());
    Lielab_domain_se.def(py::self /= int());
    Lielab_domain_se.def(py::self /= double());

    // Mic python funcs
    Lielab_domain_se.def(py::pickle([](const Lielab::domain::se& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point, obj._shape);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 2)
                throw std::runtime_error("se: Invalid state.");

            Lielab::domain::se obj;
            obj.point = t[0].cast<Lielab::domain::se::point_t>();
            obj._shape = t[1].cast<int>();
            return obj;
        }));
}
