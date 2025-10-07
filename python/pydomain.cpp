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

#include "pydomain.hpp"

namespace py = pybind11;

std::string matstr(const Eigen::VectorXd &mat);
std::string matstr(const Eigen::VectorXcd &mat);
std::string matstr(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mat);
std::string matstr(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& mat);
std::string matstr(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& mat);
std::string matstr(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& mat);

void bind_domain(py::module &m_domain)
{
    /*!
    * Bindings for Lielab::domain::cn
    */

    auto Lielab_domain_cn = py::class_<Lielab::domain::cn>(m_domain, "cn");
    Lielab_domain_cn.def_readonly_static("abelian", &Lielab::domain::cn::abelian);
    Lielab_domain_cn.def_readwrite("data", &Lielab::domain::cn::data);
    Lielab_domain_cn.def("to_string", &Lielab::domain::cn::to_string);
    Lielab_domain_cn.def(py::init<>());
    Lielab_domain_cn.def(py::init<const size_t>());
    Lielab_domain_cn.def(py::init<const Eigen::MatrixXcd&>());
    Lielab_domain_cn.def("from_shape", &Lielab::domain::cn::from_shape);
    Lielab_domain_cn.def("basis", &Lielab::domain::cn::basis);
    Lielab_domain_cn.def("project", &Lielab::domain::cn::project);
    Lielab_domain_cn.def("get_dimension", &Lielab::domain::cn::get_dimension);
    Lielab_domain_cn.def("get_shape", &Lielab::domain::cn::get_shape);
    Lielab_domain_cn.def("get_vector", &Lielab::domain::cn::get_vector);
    Lielab_domain_cn.def("get_matrix", &Lielab::domain::cn::get_matrix);
    Lielab_domain_cn.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::cn::set_vector));
    Lielab_domain_cn.def("__getitem__",
        [](const Lielab::domain::cn& self, const ptrdiff_t index)
        {
            return self[index];
        });
    Lielab_domain_cn.def("__call__", [](Lielab::domain::cn& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_cn.def("__call__", [](const Lielab::domain::cn& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_cn.def(py::self + py::self);
    Lielab_domain_cn.def(py::self += py::self);
    Lielab_domain_cn.def(py::self - py::self);
    Lielab_domain_cn.def(py::self -= py::self);
    Lielab_domain_cn.def(-py::self);
    Lielab_domain_cn.def(py::self * int());
    Lielab_domain_cn.def(py::self * double());
    // Lielab_domain_cn.def(py::self * std::complex<int>()); // TODO: Complex integers seem bugged in Eigen right now
    Lielab_domain_cn.def(py::self * std::complex<double>());
    Lielab_domain_cn.def(py::self *= int());
    Lielab_domain_cn.def(py::self *= double());
    // Lielab_domain_cn.def(py::self *= std::complex<int>()); // TODO: Complex integers seem bugged in Eigen right now
    Lielab_domain_cn.def(py::self *= std::complex<double>());
    Lielab_domain_cn.def(py::self / int());
    Lielab_domain_cn.def(py::self / double());
    // Lielab_domain_cn.def(py::self / std::complex<int>());
    Lielab_domain_cn.def(py::self / std::complex<double>());
    Lielab_domain_cn.def(py::self /= int());
    Lielab_domain_cn.def(py::self /= double());
    // Lielab_domain_cn.def(py::self /= std::complex<int>());
    Lielab_domain_cn.def(py::self /= std::complex<double>());
    Lielab_domain_cn.def(int() * py::self);
    Lielab_domain_cn.def(double() * py::self);
    // Lielab_domain_cn.def(std::complex<int>() * py::self);
    Lielab_domain_cn.def(std::complex<double>() * py::self);
    // Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    // Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::glr  & rhs) {return lhs+rhs;});
    // Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    // Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    // Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    // Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    // Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    // Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::GLR  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    // Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    // Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::glr  & rhs) {return lhs-rhs;});
    // Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    // Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    // Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    // Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    // Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    // Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});
    Lielab_domain_cn.def("__repr__", [](const Lielab::domain::cn& self)
        {
            return "<lielab.domain.cn>";
        });
    Lielab_domain_cn.def("__str__", [](const Lielab::domain::cn& self)
        {
            return matstr(self.data);
        });
    Lielab_domain_cn.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::cn::from_vector));
    Lielab_domain_cn.def_static("from_complex_vector", py::overload_cast<const Eigen::VectorXcd&>(&Lielab::domain::cn::from_complex_vector));
    Lielab_domain_cn.def("to_complex_vector", &Lielab::domain::cn::to_complex_vector);

    /*!
    * Bindings for Lielab::domain::glr
    */

    auto Lielab_domain_glr = py::class_<Lielab::domain::glr>(m_domain, "glr");
    Lielab_domain_glr.def_readonly_static("abelian", &Lielab::domain::glr::abelian);
    Lielab_domain_glr.def_readwrite("data", &Lielab::domain::glr::data);
    Lielab_domain_glr.def("to_string", &Lielab::domain::glr::to_string);
    Lielab_domain_glr.def(py::init<>());
    Lielab_domain_glr.def(py::init<const size_t>());
    Lielab_domain_glr.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_glr.def("from_shape", &Lielab::domain::glr::from_shape);
    Lielab_domain_glr.def("basis", &Lielab::domain::glr::basis);
    Lielab_domain_glr.def("project", &Lielab::domain::glr::project);
    Lielab_domain_glr.def("get_dimension", &Lielab::domain::glr::get_dimension);
    Lielab_domain_glr.def("get_shape", &Lielab::domain::glr::get_shape);
    Lielab_domain_glr.def("get_vector", &Lielab::domain::glr::get_vector);
    Lielab_domain_glr.def("get_matrix", &Lielab::domain::glr::get_matrix);
    Lielab_domain_glr.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::glr::set_vector));
    Lielab_domain_glr.def("__call__", [](const Lielab::domain::glr& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_glr.def("__call__", [](const Lielab::domain::glr& self, const ptrdiff_t index1, const ptrdiff_t index2)
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

    // Lielab_domain_gl.def("__add__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    // Lielab_domain_gl.def("__add__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::glr  & rhs) {return lhs+rhs;});
    // Lielab_domain_gl.def("__add__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    // Lielab_domain_gl.def("__add__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    // Lielab_domain_gl.def("__add__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    // Lielab_domain_gl.def("__add__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    // Lielab_domain_gl.def("__add__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    // Lielab_domain_gl.def("__add__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::GLR  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__mul__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    // Lielab_domain_gl.def("__sub__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    // Lielab_domain_gl.def("__sub__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::glr  & rhs) {return lhs-rhs;});
    // Lielab_domain_gl.def("__sub__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    // Lielab_domain_gl.def("__sub__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    // Lielab_domain_gl.def("__sub__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    // Lielab_domain_gl.def("__sub__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    // Lielab_domain_gl.def("__sub__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    // Lielab_domain_gl.def("__sub__", [](const Lielab::domain::glr  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});

    Lielab_domain_glr.def("__repr__", [](const Lielab::domain::glr& self)
        {
            return "<lielab.domain.glr>";
        });
    Lielab_domain_glr.def("__str__", [](const Lielab::domain::glr& self)
        {
            return matstr(self.data);
        });
    Lielab_domain_glr.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::glr::from_vector));
    

    /*!
    * Bindings for Lielab::domain::glc
    */

    auto Lielab_domain_glc = py::class_<Lielab::domain::glc>(m_domain, "glc");
    Lielab_domain_glc.def_readonly_static("abelian", &Lielab::domain::glc::abelian);
    Lielab_domain_glc.def_readwrite("data", &Lielab::domain::glc::data);
    Lielab_domain_glc.def("to_string", &Lielab::domain::glc::to_string);
    Lielab_domain_glc.def(py::init<>());
    Lielab_domain_glc.def(py::init<const size_t>());
    Lielab_domain_glc.def(py::init<const Eigen::MatrixXcd&>());
    Lielab_domain_glc.def("from_shape", &Lielab::domain::glc::from_shape);
    Lielab_domain_glc.def("basis", &Lielab::domain::glc::basis);
    Lielab_domain_glc.def("project", &Lielab::domain::glc::project);
    Lielab_domain_glc.def("get_dimension", &Lielab::domain::glc::get_dimension);
    Lielab_domain_glc.def("get_shape", &Lielab::domain::glc::get_shape);
    Lielab_domain_glc.def("get_vector", &Lielab::domain::glc::get_vector);
    Lielab_domain_glc.def("get_matrix", &Lielab::domain::glc::get_matrix);
    Lielab_domain_glc.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::glc::set_vector));
    Lielab_domain_glc.def("__call__", [](const Lielab::domain::glc& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_glc.def("__call__", [](const Lielab::domain::glc& self, const ptrdiff_t index1, const ptrdiff_t index2)
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
        // .def(py::self * std::complex<int>()); // TODO: Complex integers seem bugged in Eigen right now
    Lielab_domain_glc.def(py::self * std::complex<double>());
    Lielab_domain_glc.def(py::self *= int());
    Lielab_domain_glc.def(py::self *= double());
        // .def(py::self *= std::complex<int>()); // TODO: Complex integers seem bugged in Eigen right now
    Lielab_domain_glc.def(py::self *= std::complex<double>());
    Lielab_domain_glc.def(py::self / int());
    Lielab_domain_glc.def(py::self / double());
        // .def(py::self / std::complex<int>());
    Lielab_domain_glc.def(py::self / std::complex<double>());
    Lielab_domain_glc.def(py::self /= int());
    Lielab_domain_glc.def(py::self /= double());
        // .def(py::self /= std::complex<int>());
    Lielab_domain_glc.def(py::self /= std::complex<double>());
    Lielab_domain_glc.def(int() * py::self);
    Lielab_domain_glc.def(double() * py::self);
        // .def(std::complex<int>() * py::self);
    Lielab_domain_glc.def(std::complex<double>() * py::self);
    Lielab_domain_glc.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::glc::from_vector));
    Lielab_domain_glc.def_static("from_complex_vector", py::overload_cast<const Eigen::VectorXcd&>(&Lielab::domain::glc::from_complex_vector));

    // Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    // Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::glr  & rhs) {return lhs+rhs;});
    // Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    // Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    // Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    // Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    // Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    // Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::GLR  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    // Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    // Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::glr  & rhs) {return lhs-rhs;});
    // Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    // Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    // Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    // Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    // Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    // Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});

    Lielab_domain_glc.def("__repr__", [](const Lielab::domain::glc& self)
        {
            return "<lielab.domain.glc>";
        });
    Lielab_domain_glc.def("__str__", [](const Lielab::domain::glc& self)
        {
            return matstr(self.data);
        });
    

    /*!
    * Bindings for Lielab::domain::rn
    */

    auto Lielab_domain_rn = py::class_<Lielab::domain::rn>(m_domain, "rn");
    Lielab_domain_rn.def_readwrite("data", &Lielab::domain::rn::data);
    Lielab_domain_rn.def_readonly_static("abelian", &Lielab::domain::rn::abelian);
    Lielab_domain_rn.def("to_string", &Lielab::domain::rn::to_string);
    Lielab_domain_rn.def(py::init<>());
    Lielab_domain_rn.def(py::init<const size_t>());
    Lielab_domain_rn.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_rn.def("from_shape", &Lielab::domain::rn::from_shape);
    Lielab_domain_rn.def("basis", &Lielab::domain::rn::basis);
    Lielab_domain_rn.def("project", &Lielab::domain::rn::project);
    Lielab_domain_rn.def("get_dimension", &Lielab::domain::rn::get_dimension);
    Lielab_domain_rn.def("get_shape", &Lielab::domain::rn::get_shape);
    Lielab_domain_rn.def("get_vector", &Lielab::domain::rn::get_vector);
    Lielab_domain_rn.def("get_matrix", &Lielab::domain::rn::get_matrix);
    Lielab_domain_rn.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::rn::set_vector));
    Lielab_domain_rn.def("__call__", [](const Lielab::domain::rn& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_rn.def("__call__", [](const Lielab::domain::rn& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_rn.def(py::self + py::self);
    Lielab_domain_rn.def(py::self += py::self);
    Lielab_domain_rn.def(py::self - py::self);
    Lielab_domain_rn.def(py::self -= py::self);
    Lielab_domain_rn.def(-py::self);
    Lielab_domain_rn.def(py::self * int());
    Lielab_domain_rn.def(py::self * double());
    Lielab_domain_rn.def(py::self *= int());
    Lielab_domain_rn.def(py::self *= double());
    Lielab_domain_rn.def(py::self / int());
    Lielab_domain_rn.def(py::self / double());
    Lielab_domain_rn.def(py::self /= int());
    Lielab_domain_rn.def(py::self /= double());
    Lielab_domain_rn.def(int() * py::self);
    Lielab_domain_rn.def(double() * py::self);

    // Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    // Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::glr  & rhs) {return lhs+rhs;});
    // Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    // Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    // Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    // Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    // Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    // Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::GLR  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    // Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    // Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::glr  & rhs) {return lhs-rhs;});
    // Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    // Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    // Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    // Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    // Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    // Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});

    Lielab_domain_rn.def("__repr__", [](const Lielab::domain::rn& self)
        {
            return "<lielab.domain.rn>";
        });
    Lielab_domain_rn.def("__str__", [](const Lielab::domain::rn& self)
        {
            return matstr(self.data);
        });
    Lielab_domain_rn.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::rn::from_vector));
    
    /*!
    * Bindings for Lielab::domain::se
    */

    auto Lielab_domain_se = py::class_<Lielab::domain::se>(m_domain, "se");
    Lielab_domain_se.def_readwrite("data", &Lielab::domain::se::data);
    Lielab_domain_se.def_readonly_static("abelian", &Lielab::domain::se::abelian);
    Lielab_domain_se.def("to_string", &Lielab::domain::se::to_string);
    Lielab_domain_se.def(py::init<>());
    Lielab_domain_se.def(py::init<const size_t>());
    Lielab_domain_se.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_se.def("from_shape", &Lielab::domain::se::from_shape);
    Lielab_domain_se.def("basis", &Lielab::domain::se::basis);
        // .def("project", &Lielab::domain::se::project);
    Lielab_domain_se.def("get_dimension", &Lielab::domain::se::get_dimension);
    Lielab_domain_se.def("get_shape", &Lielab::domain::se::get_shape);
    Lielab_domain_se.def("get_vector", &Lielab::domain::se::get_vector);
    Lielab_domain_se.def("get_matrix", &Lielab::domain::se::get_matrix);
    Lielab_domain_se.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::se::set_vector));
    Lielab_domain_se.def("__call__", [](const Lielab::domain::se& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_se.def("__call__", [](const Lielab::domain::se& self, const ptrdiff_t index1, const ptrdiff_t index2)
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
    Lielab_domain_se.def(py::self *= int());
    Lielab_domain_se.def(py::self *= double());
    Lielab_domain_se.def(py::self / int());
    Lielab_domain_se.def(py::self / double());
    Lielab_domain_se.def(py::self /= int());
    Lielab_domain_se.def(py::self /= double());
    Lielab_domain_se.def(int() * py::self);
    Lielab_domain_se.def(double() * py::self);

    // Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    // Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::glr  & rhs) {return lhs+rhs;});
    // Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    // Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    // Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    // Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    // Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    // Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::GLR  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    // Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    // Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::glr  & rhs) {return lhs-rhs;});
    // Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    // Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    // Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    // Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    // Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    // Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});

    Lielab_domain_se.def("__repr__", [](const Lielab::domain::se& self)
        {
            return "<lielab.domain.se>";
        });
    Lielab_domain_se.def("__str__", [](const Lielab::domain::se& self)
        {
            return matstr(self.data);
        });
    Lielab_domain_se.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::se::from_vector));


    /*!
    * Bindings for Lielab::domain::so
    */

    auto Lielab_domain_so = py::class_<Lielab::domain::so>(m_domain, "so");
    Lielab_domain_so.def_readwrite("data", &Lielab::domain::so::data);
    Lielab_domain_so.def_readonly_static("abelian", &Lielab::domain::so::abelian);
    Lielab_domain_so.def("to_string", &Lielab::domain::so::to_string);
    Lielab_domain_so.def(py::init<>());
    Lielab_domain_so.def(py::init<const size_t>());
    Lielab_domain_so.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_so.def("from_shape", &Lielab::domain::so::from_shape);
    Lielab_domain_so.def("basis", &Lielab::domain::so::basis);
    Lielab_domain_so.def("project", &Lielab::domain::so::project);
    Lielab_domain_so.def("get_dimension", &Lielab::domain::so::get_dimension);
    Lielab_domain_so.def("get_shape", &Lielab::domain::so::get_shape);
    Lielab_domain_so.def("get_vector", &Lielab::domain::so::get_vector);
    Lielab_domain_so.def("get_matrix", &Lielab::domain::so::get_matrix);
    Lielab_domain_so.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::so::set_vector));
    Lielab_domain_so.def("__call__", [](const Lielab::domain::so& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_so.def("__call__", [](const Lielab::domain::so& self, const ptrdiff_t index1, const ptrdiff_t index2)
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
    Lielab_domain_so.def(py::self *= int());
    Lielab_domain_so.def(py::self *= double());
    Lielab_domain_so.def(py::self / int());
    Lielab_domain_so.def(py::self / double());
    Lielab_domain_so.def(py::self /= int());
    Lielab_domain_so.def(py::self /= double());
    Lielab_domain_so.def(int() * py::self);
    Lielab_domain_so.def(double() * py::self);

    // Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    // Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::glr  & rhs) {return lhs+rhs;});
    // Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    // Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    // Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    // Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    // Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    // Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::GLR  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    // Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    // Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::glr  & rhs) {return lhs-rhs;});
    // Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    // Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    // Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    // Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    // Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    // Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});

    Lielab_domain_so.def("__repr__", [](const Lielab::domain::so& self)
        {
            return "<lielab.domain.so>";
        });
    Lielab_domain_so.def("__str__", [](const Lielab::domain::so& self)
        {
            return matstr(self.data);
        });
    Lielab_domain_so.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::so::from_vector));
    
    /*!
    * Bindings for Lielab::domain::sp
    */
    
    auto Lielab_domain_sp = py::class_<Lielab::domain::sp>(m_domain, "sp");
    Lielab_domain_sp.def_readwrite("data", &Lielab::domain::sp::data);
    Lielab_domain_sp.def_readonly_static("abelian", &Lielab::domain::sp::abelian);
    Lielab_domain_sp.def("to_string", &Lielab::domain::sp::to_string);
    Lielab_domain_sp.def(py::init<>());
    Lielab_domain_sp.def(py::init<const size_t>());
    Lielab_domain_sp.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_sp.def("from_shape", &Lielab::domain::sp::from_shape);
    Lielab_domain_sp.def("basis", &Lielab::domain::sp::basis);
    Lielab_domain_sp.def("project", &Lielab::domain::sp::project);
    Lielab_domain_sp.def("get_dimension", &Lielab::domain::sp::get_dimension);
    Lielab_domain_sp.def("get_shape", &Lielab::domain::sp::get_shape);
    Lielab_domain_sp.def("get_vector", &Lielab::domain::sp::get_vector);
    Lielab_domain_sp.def("get_matrix", &Lielab::domain::sp::get_matrix);
    Lielab_domain_sp.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::sp::set_vector));
    Lielab_domain_sp.def("__call__", [](const Lielab::domain::sp& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_sp.def("__call__", [](const Lielab::domain::sp& self, const ptrdiff_t index1, const ptrdiff_t index2)
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
    Lielab_domain_sp.def(py::self *= int());
    Lielab_domain_sp.def(py::self *= double());
    Lielab_domain_sp.def(py::self / int());
    Lielab_domain_sp.def(py::self / double());
    Lielab_domain_sp.def(py::self /= int());
    Lielab_domain_sp.def(py::self /= double());
    Lielab_domain_sp.def(int() * py::self);
    Lielab_domain_sp.def(double() * py::self);

    // Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    // Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::glr  & rhs) {return lhs+rhs;});
    // Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    // Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    // Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    // Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    // Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    // Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::GLR  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    // Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    // Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::glr  & rhs) {return lhs-rhs;});
    // Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    // Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    // Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    // Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    // Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    // Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});

    Lielab_domain_sp.def("__repr__", [](const Lielab::domain::sp& self)
        {
            return "<lielab.domain.sp>";
        });
    Lielab_domain_sp.def("__str__", [](const Lielab::domain::sp& self)
        {
            return matstr(self.data);
        });
    Lielab_domain_sp.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::sp::from_vector));
    
    
    /*!
    * Bindings for Lielab::domain::su
    */
    
    auto Lielab_domain_su = py::class_<Lielab::domain::su>(m_domain, "su");
    Lielab_domain_su.def_readwrite("data", &Lielab::domain::su::data);
    Lielab_domain_su.def_readonly_static("abelian", &Lielab::domain::su::abelian);
    Lielab_domain_su.def("to_string", &Lielab::domain::su::to_string);
    Lielab_domain_su.def(py::init<>());
    Lielab_domain_su.def(py::init<const size_t>());
    Lielab_domain_su.def(py::init<const Eigen::MatrixXcd&>());
    Lielab_domain_su.def("from_shape", &Lielab::domain::su::from_shape);
    Lielab_domain_su.def("basis", &Lielab::domain::su::basis);
        // .def("project", &Lielab::domain::su::project); // TODO: Implement project
    Lielab_domain_su.def("get_dimension", &Lielab::domain::su::get_dimension);
    Lielab_domain_su.def("get_shape", &Lielab::domain::su::get_shape);
    Lielab_domain_su.def("get_vector", &Lielab::domain::su::get_vector);
    Lielab_domain_su.def("get_matrix", &Lielab::domain::su::get_matrix);
    Lielab_domain_su.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::su::set_vector));
    Lielab_domain_su.def("__call__", [](const Lielab::domain::su& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_su.def("__call__", [](const Lielab::domain::su& self, const ptrdiff_t index1, const ptrdiff_t index2)
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
        // .def(py::self * std::complex<int>()); // TODO: Complex integers seem bugged in Eigen right now
    Lielab_domain_su.def(py::self * std::complex<double>());
    Lielab_domain_su.def(py::self *= int());
    Lielab_domain_su.def(py::self *= double());
        // .def(py::self *= std::complex<int>()); // TODO: Complex integers seem bugged in Eigen right now
    Lielab_domain_su.def(py::self *= std::complex<double>());
    Lielab_domain_su.def(py::self / int());
    Lielab_domain_su.def(py::self / double());
        // .def(py::self / std::complex<int>());
    Lielab_domain_su.def(py::self / std::complex<double>());
    Lielab_domain_su.def(py::self /= int());
    Lielab_domain_su.def(py::self /= double());
        // .def(py::self /= std::complex<int>());
    Lielab_domain_su.def(py::self /= std::complex<double>());
    Lielab_domain_su.def(int() * py::self);
    Lielab_domain_su.def(double() * py::self);
        // .def(std::complex<int>() * py::self);
    Lielab_domain_su.def(std::complex<double>() * py::self);

    // Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    // Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::glr  & rhs) {return lhs+rhs;});
    // Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    // Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    // Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    // Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    // Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    // Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::GLR  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    // Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    // Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::glr  & rhs) {return lhs-rhs;});
    // Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    // Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    // Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    // Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    // Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    // Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});

    Lielab_domain_su.def("__repr__", [](const Lielab::domain::su& self)
        {
            return "<lielab.domain.su>";
        });
    Lielab_domain_su.def("__str__", [](const Lielab::domain::su& self)
        {
            return matstr(self.data);
        });
    Lielab_domain_su.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::su::from_vector));
    
    /*!
     * Lie Groups
     */

    auto Lielab_domain_CN = py::class_<Lielab::domain::CN>(m_domain, "CN");
    Lielab_domain_CN.def_readwrite("data", &Lielab::domain::CN::data);
    Lielab_domain_CN.def_readonly_static("abelian", &Lielab::domain::CN::abelian);
    Lielab_domain_CN.def("to_string", &Lielab::domain::CN::to_string);
    Lielab_domain_CN.def_static("from_shape", &Lielab::domain::CN::from_shape);
    Lielab_domain_CN.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::CN::from_vector));
    Lielab_domain_CN.def(py::init<>());
    Lielab_domain_CN.def(py::init<const size_t>());
    Lielab_domain_CN.def(py::init<const Eigen::MatrixXcd&>());
    // Lielab_domain_CN.def("from_shape", &Lielab::domain::CN::from_shape);
    // .def("project", &Lielab::domain::CN::project); // TODO:
    Lielab_domain_CN.def("get_dimension", &Lielab::domain::CN::get_dimension);
    Lielab_domain_CN.def("get_shape", &Lielab::domain::CN::get_shape);
    Lielab_domain_CN.def("get_matrix", &Lielab::domain::CN::get_matrix);
    Lielab_domain_CN.def("inverse", &Lielab::domain::CN::inverse);
    Lielab_domain_CN.def("serialize", &Lielab::domain::CN::serialize);
    Lielab_domain_CN.def("project", &Lielab::domain::CN::project);
    Lielab_domain_CN.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::CN::unserialize));
    Lielab_domain_CN.def("__call__", [](const Lielab::domain::CN& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_CN.def("__call__", [](const Lielab::domain::CN& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_CN.def("__getitem__",
        [](const Lielab::domain::CN& self, const ptrdiff_t index)
        {
            return self[index];
        });
    Lielab_domain_CN.def(py::self * py::self);
    Lielab_domain_CN.def(py::self *= py::self);

    // Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});

    Lielab_domain_CN.def("__repr__", [](const Lielab::domain::CN& self)
        {
            return "<lielab.domain.CN>";
        });
    Lielab_domain_CN.def("__str__", [](const Lielab::domain::CN& self)
        {
            return matstr(self.data);
        });
    Lielab_domain_CN.def_static("from_complex_vector", py::overload_cast<const Eigen::VectorXcd&>(&Lielab::domain::CN::from_complex_vector));
    Lielab_domain_CN.def("to_complex_vector", &Lielab::domain::CN::to_complex_vector);
    
    auto Lielab_domain_GLR = py::class_<Lielab::domain::GLR>(m_domain, "GLR");
    Lielab_domain_GLR.def_readwrite("data", &Lielab::domain::GLR::data);
    Lielab_domain_GLR.def_readonly_static("abelian", &Lielab::domain::GLR::abelian);
    Lielab_domain_GLR.def_static("from_shape", &Lielab::domain::GLR::from_shape);
    Lielab_domain_GLR.def("to_string", &Lielab::domain::GLR::to_string);
    Lielab_domain_GLR.def(py::init<>());
    Lielab_domain_GLR.def(py::init<const size_t>());
    Lielab_domain_GLR.def(py::init<const Eigen::MatrixXd&>());
    // Lielab_domain_GLR.def("from_shape", &Lielab::domain::GLR::from_shape);
    Lielab_domain_GLR.def("project", &Lielab::domain::GLR::project);
    Lielab_domain_GLR.def("get_dimension", &Lielab::domain::GLR::get_dimension);
    Lielab_domain_GLR.def("get_shape", &Lielab::domain::GLR::get_shape);
    Lielab_domain_GLR.def("get_matrix", &Lielab::domain::GLR::get_matrix);
    Lielab_domain_GLR.def("inverse", &Lielab::domain::GLR::inverse);
    Lielab_domain_GLR.def("serialize", &Lielab::domain::GLR::serialize);
    Lielab_domain_GLR.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::GLR::unserialize));
    Lielab_domain_GLR.def("__call__", [](const Lielab::domain::GLR& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_GLR.def(py::self * py::self);
    Lielab_domain_GLR.def(py::self *= py::self);

    // Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GLR  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GLR  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GLR  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GLR  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GLR  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GLR  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GLR  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GLR  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});

    Lielab_domain_GLR.def("__repr__", [](const Lielab::domain::GLR& self)
        {
            return "<lielab.domain.GLR>";
        });
    Lielab_domain_GLR.def("__str__", [](const Lielab::domain::GLR& self)
        {
            return matstr(self.data);
        });
    
    auto Lielab_domain_GLC = py::class_<Lielab::domain::GLC>(m_domain, "GLC");
    Lielab_domain_GLC.def_readwrite("data", &Lielab::domain::GLC::data);
    Lielab_domain_GLC.def_readonly_static("abelian", &Lielab::domain::GLC::abelian);
    Lielab_domain_GLC.def("to_string", &Lielab::domain::GLC::to_string);
    Lielab_domain_GLC.def_static("from_shape", &Lielab::domain::GLC::from_shape);
    Lielab_domain_GLC.def(py::init<>());
    Lielab_domain_GLC.def(py::init<const size_t>());
    Lielab_domain_GLC.def(py::init<const Eigen::MatrixXcd&>());
    // Lielab_domain_GLC.def("from_shape", &Lielab::domain::GLC::from_shape);
    Lielab_domain_GLC.def("project", &Lielab::domain::GLC::project);
    Lielab_domain_GLC.def("get_dimension", &Lielab::domain::GLC::get_dimension);
    Lielab_domain_GLC.def("get_shape", &Lielab::domain::GLC::get_shape);
    Lielab_domain_GLC.def("get_matrix", &Lielab::domain::GLC::get_matrix);
    Lielab_domain_GLC.def("inverse", &Lielab::domain::GLC::inverse);
    Lielab_domain_GLC.def("serialize", &Lielab::domain::GLC::serialize);
    Lielab_domain_GLC.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::GLC::unserialize));
    Lielab_domain_GLC.def("__call__", [](const Lielab::domain::GLC& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_GLC.def(py::self * py::self);
    Lielab_domain_GLC.def(py::self *= py::self);

    // Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});

    Lielab_domain_GLC.def("__repr__", [](const Lielab::domain::GLC& self)
        {
            return "<lielab.domain.GLC>";
        });
    Lielab_domain_GLC.def("__str__", [](const Lielab::domain::GLC& self)
        {
            return matstr(self.data);
        });

    auto Lielab_domain_RN = py::class_<Lielab::domain::RN>(m_domain, "RN");
    Lielab_domain_RN.def_readwrite("data", &Lielab::domain::RN::data);
    Lielab_domain_RN.def_readonly_static("abelian", &Lielab::domain::RN::abelian);
    Lielab_domain_RN.def_static("from_shape", &Lielab::domain::RN::from_shape);
    Lielab_domain_RN.def("to_string", &Lielab::domain::RN::to_string);
    Lielab_domain_RN.def(py::init<>());
    Lielab_domain_RN.def(py::init<const size_t>());
    Lielab_domain_RN.def(py::init<const Eigen::MatrixXd&>());
    // Lielab_domain_RN.def("from_shape", &Lielab::domain::RN::from_shape);
    Lielab_domain_RN.def("project", &Lielab::domain::RN::project);
    Lielab_domain_RN.def("get_dimension", &Lielab::domain::RN::get_dimension);
    Lielab_domain_RN.def("get_shape", &Lielab::domain::RN::get_shape);
    Lielab_domain_RN.def("get_matrix", &Lielab::domain::RN::get_matrix);
    Lielab_domain_RN.def("inverse", &Lielab::domain::RN::inverse);
    Lielab_domain_RN.def("serialize", &Lielab::domain::RN::serialize);
    Lielab_domain_RN.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::RN::unserialize));
    Lielab_domain_RN.def("__call__", [](const Lielab::domain::RN& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_RN.def("__call__", [](const Lielab::domain::RN& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_RN.def(py::self * py::self);
    Lielab_domain_RN.def(py::self *= py::self);

    // Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});

    Lielab_domain_RN.def("__repr__", [](const Lielab::domain::RN& self)
        {
            return "<lielab.domain.RN>";
        });
    Lielab_domain_RN.def("__str__", [](const Lielab::domain::RN& self)
        {
            return matstr(self.data);
        });
    Lielab_domain_RN.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::RN::from_vector));

    auto Lielab_domain_SE = py::class_<Lielab::domain::SE>(m_domain, "SE");
    Lielab_domain_SE.def_readwrite("data", &Lielab::domain::SE::data);
    Lielab_domain_SE.def_readonly_static("abelian", &Lielab::domain::SE::abelian);
    Lielab_domain_SE.def_static("from_shape", &Lielab::domain::SE::from_shape);
    Lielab_domain_SE.def("to_string", &Lielab::domain::SE::to_string);
    Lielab_domain_SE.def(py::init<>());
    Lielab_domain_SE.def(py::init<const size_t>());
    Lielab_domain_SE.def(py::init<const Eigen::MatrixXd&>());
    // Lielab_domain_SE.def("from_shape", &Lielab::domain::SE::from_shape);
    // .def("project", &Lielab::domain::SE::project);
    Lielab_domain_SE.def("get_dimension", &Lielab::domain::SE::get_dimension);
    Lielab_domain_SE.def("get_shape", &Lielab::domain::SE::get_shape);
    Lielab_domain_SE.def("get_matrix", &Lielab::domain::SE::get_matrix);
    Lielab_domain_SE.def("inverse", &Lielab::domain::SE::inverse);
    Lielab_domain_SE.def("serialize", &Lielab::domain::SE::serialize);
    Lielab_domain_SE.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::SE::unserialize));
    Lielab_domain_SE.def("__call__", [](const Lielab::domain::SE& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SE.def(py::self * py::self);
    Lielab_domain_SE.def(py::self *= py::self);

    // Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});

    Lielab_domain_SE.def("__repr__", [](const Lielab::domain::SE& self)
        {
            return "<lielab.domain.SE>";
        });
    Lielab_domain_SE.def("__str__", [](const Lielab::domain::SE& self)
        {
            return matstr(self.data);
        });

    auto Lielab_domain_SO = py::class_<Lielab::domain::SO>(m_domain, "SO");
    Lielab_domain_SO.def_readwrite("data", &Lielab::domain::SO::data);
    Lielab_domain_SO.def_readonly_static("abelian", &Lielab::domain::SO::abelian);
    Lielab_domain_SO.def_static("from_shape", &Lielab::domain::SO::from_shape);
    Lielab_domain_SO.def("to_string", &Lielab::domain::SO::to_string);
    Lielab_domain_SO.def(py::init<>());
    Lielab_domain_SO.def(py::init<const size_t>());
    Lielab_domain_SO.def(py::init<const Eigen::MatrixXd&>());
    // Lielab_domain_SO.def("from_shape", &Lielab::domain::SO::from_shape);
    Lielab_domain_SO.def("project", &Lielab::domain::SO::project);
    Lielab_domain_SO.def("get_dimension", &Lielab::domain::SO::get_dimension);
    Lielab_domain_SO.def("get_shape", &Lielab::domain::SO::get_shape);
    Lielab_domain_SO.def("get_matrix", &Lielab::domain::SO::get_matrix);
    Lielab_domain_SO.def("inverse", &Lielab::domain::SO::inverse);
    Lielab_domain_SO.def("serialize", &Lielab::domain::SO::serialize);
    Lielab_domain_SO.def("unserialize", py::overload_cast<const Eigen::VectorXd &>(&Lielab::domain::SO::unserialize));
    Lielab_domain_SO.def("__call__", [](const Lielab::domain::SO & self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SO.def(py::self * py::self);
    Lielab_domain_SO.def(py::self *= py::self);

    // Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});

    Lielab_domain_SO.def("__repr__", [](const Lielab::domain::SO& self)
        {
            return "<lielab.domain.SO>";
        });
    Lielab_domain_SO.def("__str__", [](const Lielab::domain::SO& self)
        {
            return matstr(self.data);
        });
    Lielab_domain_SO.def_static("from_eulerangles_body123", &Lielab::domain::SO::from_eulerangles_body123<double>, "The from_eulerangles_body123 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body231", &Lielab::domain::SO::from_eulerangles_body231<double>, "The from_eulerangles_body231 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body312", &Lielab::domain::SO::from_eulerangles_body312<double>, "The from_eulerangles_body312 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body132", &Lielab::domain::SO::from_eulerangles_body132<double>, "The from_eulerangles_body132 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body213", &Lielab::domain::SO::from_eulerangles_body213<double>, "The from_eulerangles_body213 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body321", &Lielab::domain::SO::from_eulerangles_body321<double>, "The from_eulerangles_body321 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body121", &Lielab::domain::SO::from_eulerangles_body121<double>, "The from_eulerangles_body121 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body131", &Lielab::domain::SO::from_eulerangles_body131<double>, "The from_eulerangles_body131 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body212", &Lielab::domain::SO::from_eulerangles_body212<double>, "The from_eulerangles_body212 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body232", &Lielab::domain::SO::from_eulerangles_body232<double>, "The from_eulerangles_body232 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body313", &Lielab::domain::SO::from_eulerangles_body313<double>, "The from_eulerangles_body313 function.");
    Lielab_domain_SO.def_static("from_eulerangles_body323", &Lielab::domain::SO::from_eulerangles_body323<double>, "The from_eulerangles_body323 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space123", &Lielab::domain::SO::from_eulerangles_space123<double>, "The from_eulerangles_space123 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space231", &Lielab::domain::SO::from_eulerangles_space231<double>, "The from_eulerangles_space231 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space312", &Lielab::domain::SO::from_eulerangles_space312<double>, "The from_eulerangles_space312 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space132", &Lielab::domain::SO::from_eulerangles_space132<double>, "The from_eulerangles_space132 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space213", &Lielab::domain::SO::from_eulerangles_space213<double>, "The from_eulerangles_space213 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space321", &Lielab::domain::SO::from_eulerangles_space321<double>, "The from_eulerangles_space321 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space121", &Lielab::domain::SO::from_eulerangles_space121<double>, "The from_eulerangles_space121 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space131", &Lielab::domain::SO::from_eulerangles_space131<double>, "The from_eulerangles_space131 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space212", &Lielab::domain::SO::from_eulerangles_space212<double>, "The from_eulerangles_space212 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space232", &Lielab::domain::SO::from_eulerangles_space232<double>, "The from_eulerangles_space232 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space313", &Lielab::domain::SO::from_eulerangles_space313<double>, "The from_eulerangles_space313 function.");
    Lielab_domain_SO.def_static("from_eulerangles_space323", &Lielab::domain::SO::from_eulerangles_space323<double>, "The from_eulerangles_space323 function.");
    Lielab_domain_SO.def_static("from_quaternion", &Lielab::domain::SO::from_quaternion<double>, "The from_quaternion function.");
    Lielab_domain_SO.def_static("from_rodriguesvector", &Lielab::domain::SO::from_rodriguesvector<double>, "The from_rodriguesvector function.");
    Lielab_domain_SO.def_static("from_SU2", &Lielab::domain::SO::from_SU2);
    Lielab_domain_SO.def("to_quaternion", &Lielab::domain::SO::to_quaternion<double>);
    Lielab_domain_SO.def("to_gibbs", &Lielab::domain::SO::to_gibbs<double>);
    Lielab_domain_SO.def("to_eulerangles_body123", &Lielab::domain::SO::to_eulerangles_body123<double>, "The to_eulerangles_body123 function.");
    Lielab_domain_SO.def("to_eulerangles_body231", &Lielab::domain::SO::to_eulerangles_body231<double>, "The to_eulerangles_body231 function.");
    Lielab_domain_SO.def("to_eulerangles_body312", &Lielab::domain::SO::to_eulerangles_body312<double>, "The to_eulerangles_body312 function.");
    Lielab_domain_SO.def("to_eulerangles_body132", &Lielab::domain::SO::to_eulerangles_body132<double>, "The to_eulerangles_body132 function.");
    Lielab_domain_SO.def("to_eulerangles_body213", &Lielab::domain::SO::to_eulerangles_body213<double>, "The to_eulerangles_body213 function.");
    Lielab_domain_SO.def("to_eulerangles_body321", &Lielab::domain::SO::to_eulerangles_body321<double>, "The to_eulerangles_body321 function.");
    Lielab_domain_SO.def("to_eulerangles_body121", &Lielab::domain::SO::to_eulerangles_body121<double>, "The to_eulerangles_body121 function.");
    Lielab_domain_SO.def("to_eulerangles_body131", &Lielab::domain::SO::to_eulerangles_body131<double>, "The to_eulerangles_body131 function.");
    Lielab_domain_SO.def("to_eulerangles_body212", &Lielab::domain::SO::to_eulerangles_body212<double>, "The to_eulerangles_body212 function.");
    Lielab_domain_SO.def("to_eulerangles_body232", &Lielab::domain::SO::to_eulerangles_body232<double>, "The to_eulerangles_body232 function.");
    Lielab_domain_SO.def("to_eulerangles_body313", &Lielab::domain::SO::to_eulerangles_body313<double>, "The to_eulerangles_body313 function.");
    Lielab_domain_SO.def("to_eulerangles_body323", &Lielab::domain::SO::to_eulerangles_body323<double>, "The to_eulerangles_body323 function.");
    Lielab_domain_SO.def("to_eulerangles_space123", &Lielab::domain::SO::to_eulerangles_space123<double>, "The to_eulerangles_space123 function.");
    Lielab_domain_SO.def("to_eulerangles_space231", &Lielab::domain::SO::to_eulerangles_space231<double>, "The to_eulerangles_space231 function.");
    Lielab_domain_SO.def("to_eulerangles_space312", &Lielab::domain::SO::to_eulerangles_space312<double>, "The to_eulerangles_space312 function.");
    Lielab_domain_SO.def("to_eulerangles_space132", &Lielab::domain::SO::to_eulerangles_space132<double>, "The to_eulerangles_space132 function.");
    Lielab_domain_SO.def("to_eulerangles_space213", &Lielab::domain::SO::to_eulerangles_space213<double>, "The to_eulerangles_space213 function.");
    Lielab_domain_SO.def("to_eulerangles_space321", &Lielab::domain::SO::to_eulerangles_space321<double>, "The to_eulerangles_space321 function.");
    Lielab_domain_SO.def("to_eulerangles_space121", &Lielab::domain::SO::to_eulerangles_space121<double>, "The to_eulerangles_space121 function.");
    Lielab_domain_SO.def("to_eulerangles_space131", &Lielab::domain::SO::to_eulerangles_space131<double>, "The to_eulerangles_space131 function.");
    Lielab_domain_SO.def("to_eulerangles_space212", &Lielab::domain::SO::to_eulerangles_space212<double>, "The to_eulerangles_space212 function.");
    Lielab_domain_SO.def("to_eulerangles_space232", &Lielab::domain::SO::to_eulerangles_space232<double>, "The to_eulerangles_space232 function.");
    Lielab_domain_SO.def("to_eulerangles_space313", &Lielab::domain::SO::to_eulerangles_space313<double>, "The to_eulerangles_space313 function.");
    Lielab_domain_SO.def("to_eulerangles_space323", &Lielab::domain::SO::to_eulerangles_space323<double>, "The to_eulerangles_space323 function.");
    

    auto Lielab_domain_SP = py::class_<Lielab::domain::SP>(m_domain, "SP");
    Lielab_domain_SP.def_readwrite("data", &Lielab::domain::SP::data);
    Lielab_domain_SP.def_readonly_static("abelian", &Lielab::domain::SP::abelian);
    Lielab_domain_SP.def_static("from_shape", Lielab::domain::SP::from_shape);
    Lielab_domain_SP.def("to_string", &Lielab::domain::SP::to_string);
    Lielab_domain_SP.def(py::init<>());
    Lielab_domain_SP.def(py::init<const size_t>());
    Lielab_domain_SP.def(py::init<const Eigen::MatrixXd&>());
    // Lielab_domain_SP.def("from_shape", &Lielab::domain::SP::from_shape);
    // .def("project", &Lielab::domain::SP::project); // TODO:
    Lielab_domain_SP.def("get_dimension", &Lielab::domain::SP::get_dimension);
    Lielab_domain_SP.def("get_shape", &Lielab::domain::SP::get_shape);
    Lielab_domain_SP.def("get_matrix", &Lielab::domain::SP::get_matrix);
    Lielab_domain_SP.def("inverse", &Lielab::domain::SP::inverse);
    Lielab_domain_SP.def("serialize", &Lielab::domain::SP::serialize);
    Lielab_domain_SP.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::SP::unserialize));
    Lielab_domain_SP.def("__call__", [](const Lielab::domain::SP& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SP.def(py::self * py::self);
    Lielab_domain_SP.def(py::self *= py::self);

    // Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});

    Lielab_domain_SP.def("__repr__", [](const Lielab::domain::SP& self)
        {
            return "<lielab.domain.SP>";
        });
    Lielab_domain_SP.def("__str__", [](const Lielab::domain::SP& self)
        {
            return matstr(self.data);
        });

    auto Lielab_domain_SU = py::class_<Lielab::domain::SU>(m_domain, "SU");
    Lielab_domain_SU.def_readwrite("data", &Lielab::domain::SU::data);
    Lielab_domain_SU.def_readonly_static("abelian", &Lielab::domain::SU::abelian);
    Lielab_domain_SU.def_static("from_shape", Lielab::domain::SU::from_shape);
    Lielab_domain_SU.def("to_string", &Lielab::domain::SU::to_string);
    Lielab_domain_SU.def(py::init<>());
    Lielab_domain_SU.def(py::init<const size_t>());
    Lielab_domain_SU.def(py::init<const Eigen::MatrixXcd&>());
    // Lielab_domain_SU.def("from_shape", &Lielab::domain::SU::from_shape);
    // .def("project", &Lielab::domain::SU::project); // TODO:
    Lielab_domain_SU.def("get_dimension", &Lielab::domain::SU::get_dimension);
    Lielab_domain_SU.def("get_shape", &Lielab::domain::SU::get_shape);
    Lielab_domain_SU.def("get_matrix", &Lielab::domain::SU::get_matrix);
    Lielab_domain_SU.def("inverse", &Lielab::domain::SU::inverse);
    Lielab_domain_SU.def("serialize", &Lielab::domain::SU::serialize);
    Lielab_domain_SU.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::SU::unserialize));
    Lielab_domain_SU.def("__call__", [](const Lielab::domain::SU& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SU.def(py::self * py::self);
    Lielab_domain_SU.def(py::self *= py::self);

    // Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    // Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::glr  & rhs) {return lhs*rhs;});
    // Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    // Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    // Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    // Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    // Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    // Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});

    Lielab_domain_SU.def("__repr__", [](const Lielab::domain::SU & self)
        {
            return "<lielab.domain.SU>";
        });
    Lielab_domain_SU.def("__str__", [](const Lielab::domain::SU & self)
        {
            return matstr(self.data);
        });
    Lielab_domain_SU.def_static("from_quaternion", &Lielab::domain::SU::from_quaternion<double>);
    Lielab_domain_SU.def_static("from_SO3", &Lielab::domain::SU::from_SO3);
    Lielab_domain_SU.def("to_quaternion", &Lielab::domain::SU::to_quaternion);

    auto Lielab_domain_CompositeAlgebra = py::class_<Lielab::domain::CompositeAlgebra>(m_domain, "CompositeAlgebra");
    Lielab_domain_CompositeAlgebra.def("to_string", &Lielab::domain::CompositeAlgebra::to_string);
    Lielab_domain_CompositeAlgebra.def(py::init());
    Lielab_domain_CompositeAlgebra.def(py::init<const size_t>());
    Lielab_domain_CompositeAlgebra.def(py::init<const std::vector<Lielab::domain::CompositeAlgebra::TYPES>&>());
    Lielab_domain_CompositeAlgebra.def("from_shape", &Lielab::domain::CompositeAlgebra::from_shape);
    Lielab_domain_CompositeAlgebra.def("basis", &Lielab::domain::CompositeAlgebra::basis);
    Lielab_domain_CompositeAlgebra.def_readwrite("space", &Lielab::domain::CompositeAlgebra::space);
    Lielab_domain_CompositeAlgebra.def("get_dimension", &Lielab::domain::CompositeAlgebra::get_dimension);
    Lielab_domain_CompositeAlgebra.def("get_shape", &Lielab::domain::CompositeAlgebra::get_shape);
    Lielab_domain_CompositeAlgebra.def("get_shapes", &Lielab::domain::CompositeAlgebra::get_shapes);
    Lielab_domain_CompositeAlgebra.def("get_dimensions", &Lielab::domain::CompositeAlgebra::get_dimensions);
    Lielab_domain_CompositeAlgebra.def("get_vectors", &Lielab::domain::CompositeAlgebra::get_vectors);
    Lielab_domain_CompositeAlgebra.def("get_shapes", &Lielab::domain::CompositeAlgebra::get_shapes);
    Lielab_domain_CompositeAlgebra.def("get_vector", &Lielab::domain::CompositeAlgebra::get_vector);
    Lielab_domain_CompositeAlgebra.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::CompositeAlgebra::set_vector));
    Lielab_domain_CompositeAlgebra.def("get_matrix", &Lielab::domain::CompositeAlgebra::get_matrix);
    Lielab_domain_CompositeAlgebra.def(py::self + py::self);
    Lielab_domain_CompositeAlgebra.def(py::self += py::self);
    Lielab_domain_CompositeAlgebra.def(py::self - py::self);
    Lielab_domain_CompositeAlgebra.def(py::self -= py::self);
    Lielab_domain_CompositeAlgebra.def(-py::self);
    Lielab_domain_CompositeAlgebra.def(py::self * int());
    Lielab_domain_CompositeAlgebra.def(py::self * double());
    // Lielab_domain_CompositeAlgebra.def(py::self * py::self);
    Lielab_domain_CompositeAlgebra.def(py::self *= int());
    Lielab_domain_CompositeAlgebra.def(py::self *= double());
    // Lielab_domain_CompositeAlgebra.def(py::self *= py::self);
    Lielab_domain_CompositeAlgebra.def(py::self / int());
    Lielab_domain_CompositeAlgebra.def(py::self / double());
    Lielab_domain_CompositeAlgebra.def(py::self /= int());
    Lielab_domain_CompositeAlgebra.def(py::self /= double());
    Lielab_domain_CompositeAlgebra.def(int() * py::self);
    Lielab_domain_CompositeAlgebra.def(double() * py::self);
    Lielab_domain_CompositeAlgebra.def("__call__",
        [](const Lielab::domain::CompositeAlgebra& self, const ptrdiff_t index)
        {
            return self(index);
        });
    Lielab_domain_CompositeAlgebra.def("__call__",
        [](const Lielab::domain::CompositeAlgebra& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_CompositeAlgebra.def("__getitem__",
        [](const Lielab::domain::CompositeAlgebra& self, const ptrdiff_t index)
        {
            return self[index];
        });
    Lielab_domain_CompositeAlgebra.def("__repr__",
        [](const Lielab::domain::CompositeAlgebra& self)
        {
            return "<lielab.domain.CompositeAlgebra>";
        });
    Lielab_domain_CompositeAlgebra.def("__str__",
        [](const Lielab::domain::CompositeAlgebra& self)
        {
            return self.to_string();
        });
    

    auto Lielab_domain_CompositeGroup = py::class_<Lielab::domain::CompositeGroup>(m_domain, "CompositeGroup");
    Lielab_domain_CompositeGroup.def("to_string", &Lielab::domain::CompositeGroup::to_string);
    Lielab_domain_CompositeGroup.def(py::init());
    Lielab_domain_CompositeGroup.def(py::init<const size_t>());
    Lielab_domain_CompositeGroup.def(py::init<const std::vector<Lielab::domain::CompositeGroup::TYPES>&>());
    Lielab_domain_CompositeGroup.def("from_shape", &Lielab::domain::CompositeGroup::from_shape);
    Lielab_domain_CompositeGroup.def_readwrite("space", &Lielab::domain::CompositeGroup::space);
    Lielab_domain_CompositeGroup.def("get_dimension", &Lielab::domain::CompositeGroup::get_dimension);
    Lielab_domain_CompositeGroup.def("get_dimensions", &Lielab::domain::CompositeGroup::get_dimensions);
    Lielab_domain_CompositeGroup.def("get_shape", &Lielab::domain::CompositeGroup::get_shape);
    Lielab_domain_CompositeGroup.def("get_shapes", &Lielab::domain::CompositeGroup::get_shapes);
    Lielab_domain_CompositeGroup.def("get_size", &Lielab::domain::CompositeGroup::get_size);
    Lielab_domain_CompositeGroup.def("get_sizes", &Lielab::domain::CompositeGroup::get_sizes);
    Lielab_domain_CompositeGroup.def("serialize", &Lielab::domain::CompositeGroup::serialize);
    Lielab_domain_CompositeGroup.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::CompositeGroup::unserialize));
    Lielab_domain_CompositeGroup.def("get_matrix", &Lielab::domain::CompositeGroup::get_matrix);
    Lielab_domain_CompositeGroup.def(py::self * py::self);
    Lielab_domain_CompositeGroup.def(py::self *= py::self);
    Lielab_domain_CompositeGroup.def("inverse", &Lielab::domain::CompositeGroup::inverse);
    Lielab_domain_CompositeGroup.def("__call__",
        [](const Lielab::domain::CompositeGroup& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_CompositeGroup.def("__getitem__",
        [](const Lielab::domain::CompositeGroup& self, const ptrdiff_t index)
        {
            return self[index];
        });
    Lielab_domain_CompositeGroup.def("__repr__",
        [](const Lielab::domain::CompositeGroup & self)
        {
            return "<lielab.domain.CompositeGroup>";
        });
    Lielab_domain_CompositeGroup.def("__str__",
        [](const Lielab::domain::CompositeGroup & self)
        {
            return self.to_string();
        });
    
    auto Lielab_domain_GrR = py::class_<Lielab::domain::GrR>(m_domain, "GrR");
    Lielab_domain_GrR.def_readwrite("data", &Lielab::domain::GrR::data);
    Lielab_domain_GrR.def("to_string", &Lielab::domain::GrR::to_string);
    Lielab_domain_GrR.def(py::init());
    Lielab_domain_GrR.def(py::init<const size_t, const size_t>());
    Lielab_domain_GrR.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_GrR.def("get_dimension", &Lielab::domain::GrR::get_dimension);
    Lielab_domain_GrR.def("get_size", &Lielab::domain::GrR::get_size);
    Lielab_domain_GrR.def("serialize", &Lielab::domain::GrR::serialize);
    Lielab_domain_GrR.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::GrR::unserialize));
    Lielab_domain_GrR.def("get_matrix", &Lielab::domain::GrR::get_matrix);
    Lielab_domain_GrR.def("project_onto", &Lielab::domain::GrR::project_onto);
    Lielab_domain_GrR.def("__call__",
        [](const Lielab::domain::GrR& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    // Lielab_domain_GrR.def("__getitem__",
    //     [](const Lielab::domain::GrR& self, const ptrdiff_t index)
    //     {
    //         return self[index];
    //     });
    Lielab_domain_GrR.def("__repr__",
        [](const Lielab::domain::GrR& self)
        {
            return "<lielab.domain.GrR>";
        });
    Lielab_domain_GrR.def("__str__",
        [](const Lielab::domain::GrR& self)
        {
            return self.to_string();
        });
    
    auto Lielab_domain_CompositeManifold = py::class_<Lielab::domain::CompositeManifold>(m_domain, "CompositeManifold");
    Lielab_domain_CompositeManifold.def("to_string", &Lielab::domain::CompositeManifold::to_string);
    Lielab_domain_CompositeManifold.def(py::init());
    Lielab_domain_CompositeManifold.def(py::init<const size_t>());
    Lielab_domain_CompositeManifold.def(py::init<const std::vector<Lielab::domain::CompositeManifold::TYPES>&>());
    Lielab_domain_CompositeManifold.def("from_shape", &Lielab::domain::CompositeManifold::from_shape);
    Lielab_domain_CompositeManifold.def_readwrite("space", &Lielab::domain::CompositeManifold::space);
    Lielab_domain_CompositeManifold.def("get_dimension", &Lielab::domain::CompositeManifold::get_dimension);
    Lielab_domain_CompositeManifold.def("get_dimensions", &Lielab::domain::CompositeManifold::get_dimensions);
    Lielab_domain_CompositeManifold.def("get_shape", &Lielab::domain::CompositeManifold::get_shape);
    Lielab_domain_CompositeManifold.def("get_shapes", &Lielab::domain::CompositeManifold::get_shapes);
    Lielab_domain_CompositeManifold.def("get_size", &Lielab::domain::CompositeManifold::get_size);
    Lielab_domain_CompositeManifold.def("get_sizes", &Lielab::domain::CompositeManifold::get_sizes);
    Lielab_domain_CompositeManifold.def("serialize", &Lielab::domain::CompositeManifold::serialize);
    Lielab_domain_CompositeManifold.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::CompositeManifold::unserialize));
    Lielab_domain_CompositeManifold.def("get_matrix", &Lielab::domain::CompositeManifold::get_matrix);
    Lielab_domain_CompositeManifold.def("__call__",
        [](const Lielab::domain::CompositeManifold& self, const ptrdiff_t index1, const ptrdiff_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_CompositeManifold.def("__getitem__",
        [](const Lielab::domain::CompositeManifold& self, const ptrdiff_t index)
        {
            return self[index];
        });
    Lielab_domain_CompositeManifold.def("__repr__",
        [](const Lielab::domain::CompositeManifold& self)
        {
            return "<lielab.domain.CompositeManifold>";
        });
    Lielab_domain_CompositeManifold.def("__str__",
        [](const Lielab::domain::CompositeManifold& self)
        {
            return self.to_string();
        });
}
