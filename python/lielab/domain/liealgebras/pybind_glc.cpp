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

#include "pybind_cn.hpp"

namespace py = pybind11;

void bind_glc(py::module& m_domain)
{
    auto Lielab_domain_glc = py::class_<Lielab::domain::glc>(m_domain, "glc");
    Lielab_domain_glc.def_readwrite("point", &Lielab::domain::glc::point);
    Lielab_domain_glc.def(py::init<>());
    Lielab_domain_glc.def(py::init<const Eigen::MatrixXcd&>());
    Lielab_domain_glc.def_static("basis", &Lielab::domain::glc::basis);
    Lielab_domain_glc.def_static("zero", &Lielab::domain::glc::zero);
    Lielab_domain_glc.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::glc::from_vector));
    Lielab_domain_glc.def_static("project", &Lielab::domain::glc::project);
    Lielab_domain_glc.def(py::init<const int>());
    Lielab_domain_glc.def_static("from_complex_vector", py::overload_cast<const Eigen::VectorXcd&>(&Lielab::domain::glc::from_complex_vector));
    Lielab_domain_glc.def("to_string", &Lielab::domain::glc::to_string);
    Lielab_domain_glc.def("__repr__", [](const Lielab::domain::glc& self)
        {
            return "<lielab.domain.glc>";
        });
    Lielab_domain_glc.def("__str__", [](const Lielab::domain::glc& self)
        {
            return "<lielab.domain.glc>";
        });
    Lielab_domain_glc.def("get_dimension", &Lielab::domain::glc::get_dimension);
    Lielab_domain_glc.def("get_size", &Lielab::domain::glc::get_size);
    Lielab_domain_glc.def("is_abelian", &Lielab::domain::glc::is_abelian);
    Lielab_domain_glc.def("get_shape", &Lielab::domain::glc::get_shape);
    Lielab_domain_glc.def("get_point", &Lielab::domain::glc::get_point);
    Lielab_domain_glc.def("serialize", &Lielab::domain::glc::serialize);
    Lielab_domain_glc.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::glc::unserialize));
    Lielab_domain_glc.def("get_matrix", &Lielab::domain::glc::get_matrix);
    Lielab_domain_glc.def("get_vector", &Lielab::domain::glc::get_vector);
    Lielab_domain_glc.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::glc::set_vector));
    Lielab_domain_glc.def("__call__", [](const Lielab::domain::glc& self, const int index)
        {
            return self(index);
        });
    Lielab_domain_glc.def("__call__", [](const Lielab::domain::glc& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    
    Lielab_domain_glc.def(py::self + py::self);
    Lielab_domain_glc.def(py::self += py::self);
    Lielab_domain_glc.def(py::self - py::self);
    Lielab_domain_glc.def(py::self -= py::self);
    Lielab_domain_glc.def(-py::self);
    Lielab_domain_glc.def(py::self * int());
    Lielab_domain_glc.def(py::self * double());
    Lielab_domain_glc.def(int() * py::self);
    Lielab_domain_glc.def(double() * py::self);
    Lielab_domain_glc.def(py::self *= int());
    Lielab_domain_glc.def(py::self *= double());
    Lielab_domain_glc.def(py::self / int());
    Lielab_domain_glc.def(py::self / double());
    Lielab_domain_glc.def(py::self /= int());
    Lielab_domain_glc.def(py::self /= double());
    Lielab_domain_glc.def(py::self * std::complex<double>());
    Lielab_domain_glc.def(std::complex<double>() * py::self);
    Lielab_domain_glc.def(py::self *= std::complex<double>());
    Lielab_domain_glc.def(py::self / std::complex<double>());
    Lielab_domain_glc.def(py::self /= std::complex<double>());
    
    // Other misc python
    Lielab_domain_glc.def(py::pickle([](const Lielab::domain::glc& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 1)
                throw std::runtime_error("glc: Invalid state.");

            Lielab::domain::glc obj;
            obj.point = t[0].cast<Lielab::domain::glc::point_t>();
            return obj;
        }));
}
