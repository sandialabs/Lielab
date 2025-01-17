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

namespace py = pybind11;

std::string matstr(const Eigen::VectorXd &mat){
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

std::string matstr(const Eigen::VectorXcd &mat){
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

std::string matstr(const Eigen::MatrixXd &mat){
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

std::string matstr(const Eigen::MatrixXcd &mat){
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

PYBIND11_MODULE(cppLielab, m) {
    /*!
    * Metadata.
    */

    m.doc() = "Lielab Python plugin";
    m.attr("__author__") = Lielab::AUTHOR;
    m.attr("__contact__") = Lielab::CONTACT;
    m.attr("__location__") = Lielab::LOCATION;
    m.attr("__version__") = Lielab::VERSION;

    /*!
    * Begin content for the "lielab" module.
    */

    py::enum_<Lielab::ALGO_STATUS>(m, "ALGO_STATUS")
        .value("OK", Lielab::ALGO_STATUS::OK)
        .value("MAXITER", Lielab::ALGO_STATUS::MAXITER)
        .value("FINISHED", Lielab::ALGO_STATUS::FINISHED);

    /*!
    * Begin content for the "domain" submodule.
    */
    py::module m_domain = m.def_submodule("domain", "The domain submodule.");

    /*!
    * Bindings for Lielab::domain::cn
    */

    auto Lielab_domain_cn = py::class_<Lielab::domain::cn>(m_domain, "cn");
    Lielab_domain_cn.def_readwrite("_data", &Lielab::domain::cn::_data);
    Lielab_domain_cn.def_readonly_static("abelian", &Lielab::domain::cn::abelian);
    Lielab_domain_cn.def(py::init<>());
    Lielab_domain_cn.def(py::init<const size_t>());
    Lielab_domain_cn.def(py::init<const Eigen::MatrixXcd &>());
    Lielab_domain_cn.def("basis", &Lielab::domain::cn::basis);
    Lielab_domain_cn.def("project", &Lielab::domain::cn::project);
    Lielab_domain_cn.def("get_dimension", &Lielab::domain::cn::get_dimension);
    Lielab_domain_cn.def("get_vector", &Lielab::domain::cn::get_vector);
    Lielab_domain_cn.def("get_matrix", &Lielab::domain::cn::get_matrix);
    Lielab_domain_cn.def("set_vector", &Lielab::domain::cn::set_vector);
    Lielab_domain_cn.def("__call__", [](Lielab::domain::cn & self, const size_t index)
        {
            return self(index);
        });
    Lielab_domain_cn.def("__call__", [](const Lielab::domain::cn & self, const size_t index1, const size_t index2)
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
        // .def(py::self * std::complex<int>()); // TODO: Complex integers seem bugged in Eigen right now
    Lielab_domain_cn.def(py::self * std::complex<double>());
    Lielab_domain_cn.def(py::self *= int());
    Lielab_domain_cn.def(py::self *= double());
        // .def(py::self *= std::complex<int>()); // TODO: Complex integers seem bugged in Eigen right now
    Lielab_domain_cn.def(py::self *= std::complex<double>());
    Lielab_domain_cn.def(py::self *= py::self);
    Lielab_domain_cn.def(py::self / int());
    Lielab_domain_cn.def(py::self / double());
        // .def(py::self / std::complex<int>());
    Lielab_domain_cn.def(py::self / std::complex<double>());
    Lielab_domain_cn.def(py::self /= int());
    Lielab_domain_cn.def(py::self /= double());
        // .def(py::self /= std::complex<int>());
    Lielab_domain_cn.def(py::self /= std::complex<double>());
    Lielab_domain_cn.def(int() * py::self);
    Lielab_domain_cn.def(double() * py::self);
        // .def(std::complex<int>() * py::self);
    Lielab_domain_cn.def(std::complex<double>() * py::self);
    Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::gl  & rhs) {return lhs+rhs;});
    Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    Lielab_domain_cn.def("__add__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::GL  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__mul__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::gl  & rhs) {return lhs-rhs;});
    Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    Lielab_domain_cn.def("__sub__", [](const Lielab::domain::cn  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});
    Lielab_domain_cn.def("__repr__", [](const Lielab::domain::cn & self)
        {
            return "<lielab.domain.cn>";
        });
    Lielab_domain_cn.def("__str__", [](const Lielab::domain::cn & self)
        {
            return matstr(self._data);
        });

    /*!
    * Bindings for Lielab::domain::gl
    */

    auto Lielab_domain_gl = py::class_<Lielab::domain::gl>(m_domain, "gl");
    Lielab_domain_gl.def_readwrite("_data", &Lielab::domain::gl::_data);
    Lielab_domain_gl.def_readonly_static("abelian", &Lielab::domain::gl::abelian);
    Lielab_domain_gl.def(py::init<>());
    Lielab_domain_gl.def(py::init<const size_t>());
    Lielab_domain_gl.def(py::init<const Eigen::MatrixXd &>());
    Lielab_domain_gl.def("basis", &Lielab::domain::gl::basis);
    Lielab_domain_gl.def("project", &Lielab::domain::gl::project);
    Lielab_domain_gl.def("get_dimension", &Lielab::domain::gl::get_dimension);
    Lielab_domain_gl.def("get_vector", &Lielab::domain::gl::get_vector);
    Lielab_domain_gl.def("get_matrix", &Lielab::domain::gl::get_matrix);
    Lielab_domain_gl.def("set_vector", &Lielab::domain::gl::set_vector);
    Lielab_domain_gl.def("__call__", [](Lielab::domain::gl & self, const size_t index)
        {
            return self(index);
        });
    Lielab_domain_gl.def("__call__", [](const Lielab::domain::gl & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_gl.def(py::self + py::self);
    Lielab_domain_gl.def(py::self += py::self);
    Lielab_domain_gl.def(py::self - py::self);
    Lielab_domain_gl.def(py::self -= py::self);
    Lielab_domain_gl.def(-py::self);
    Lielab_domain_gl.def(py::self * int());
    Lielab_domain_gl.def(py::self * double());
    Lielab_domain_gl.def(py::self *= int());
    Lielab_domain_gl.def(py::self *= double());
    Lielab_domain_gl.def(py::self *= py::self);
    Lielab_domain_gl.def(py::self / int());
    Lielab_domain_gl.def(py::self / double());
    Lielab_domain_gl.def(py::self /= int());
    Lielab_domain_gl.def(py::self /= double());
    Lielab_domain_gl.def(int() * py::self);
    Lielab_domain_gl.def(double() * py::self);
    Lielab_domain_gl.def("__add__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    Lielab_domain_gl.def("__add__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::gl  & rhs) {return lhs+rhs;});
    Lielab_domain_gl.def("__add__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    Lielab_domain_gl.def("__add__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    Lielab_domain_gl.def("__add__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    Lielab_domain_gl.def("__add__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    Lielab_domain_gl.def("__add__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    Lielab_domain_gl.def("__add__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::GL  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__mul__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    Lielab_domain_gl.def("__sub__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    Lielab_domain_gl.def("__sub__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::gl  & rhs) {return lhs-rhs;});
    Lielab_domain_gl.def("__sub__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    Lielab_domain_gl.def("__sub__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    Lielab_domain_gl.def("__sub__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    Lielab_domain_gl.def("__sub__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    Lielab_domain_gl.def("__sub__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    Lielab_domain_gl.def("__sub__", [](const Lielab::domain::gl  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});
    Lielab_domain_gl.def("__repr__", [](const Lielab::domain::gl & self)
        {
            return "<lielab.domain.gl>";
        });
    Lielab_domain_gl.def("__str__", [](const Lielab::domain::gl & self)
        {
            return matstr(self._data);
        });
    

    /*!
    * Bindings for Lielab::domain::glc
    */

    auto Lielab_domain_glc = py::class_<Lielab::domain::glc>(m_domain, "glc");
    Lielab_domain_glc.def_readwrite("_data", &Lielab::domain::glc::_data);
    Lielab_domain_glc.def_readonly_static("abelian", &Lielab::domain::glc::abelian);
    Lielab_domain_glc.def(py::init<>());
    Lielab_domain_glc.def(py::init<const size_t>());
    Lielab_domain_glc.def(py::init<const Eigen::MatrixXcd &>());
    Lielab_domain_glc.def("basis", &Lielab::domain::glc::basis);
    Lielab_domain_glc.def("project", &Lielab::domain::glc::project);
    Lielab_domain_glc.def("get_dimension", &Lielab::domain::glc::get_dimension);
    Lielab_domain_glc.def("get_vector", &Lielab::domain::glc::get_vector);
    Lielab_domain_glc.def("get_matrix", &Lielab::domain::glc::get_matrix);
    Lielab_domain_glc.def("set_vector", &Lielab::domain::glc::set_vector);
    Lielab_domain_glc.def("__call__", [](Lielab::domain::glc & self, const size_t index)
        {
            return self(index);
        });
    Lielab_domain_glc.def("__call__", [](const Lielab::domain::glc & self, const size_t index1, const size_t index2)
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
    Lielab_domain_glc.def(py::self *= py::self);
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
    Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::gl  & rhs) {return lhs+rhs;});
    Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    Lielab_domain_glc.def("__add__", [](const Lielab::domain::glc & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::GL  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__mul__", [](const Lielab::domain::glc & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::gl  & rhs) {return lhs-rhs;});
    Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    Lielab_domain_glc.def("__sub__", [](const Lielab::domain::glc & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});
    Lielab_domain_glc.def("__repr__", [](const Lielab::domain::glc & self)
        {
            return "<lielab.domain.glc>";
        });
    Lielab_domain_glc.def("__str__", [](const Lielab::domain::glc & self)
        {
            return matstr(self._data);
        });
    

    /*!
    * Bindings for Lielab::domain::rn
    */

    auto Lielab_domain_rn = py::class_<Lielab::domain::rn>(m_domain, "rn");
    Lielab_domain_rn.def_readwrite("_data", &Lielab::domain::rn::_data);
    Lielab_domain_rn.def_readwrite("shape", &Lielab::domain::rn::shape);
    Lielab_domain_rn.def_readonly_static("abelian", &Lielab::domain::rn::abelian);
    Lielab_domain_rn.def(py::init<>());
    Lielab_domain_rn.def(py::init<const size_t>());
    Lielab_domain_rn.def(py::init<const Eigen::MatrixXd &>());
    Lielab_domain_rn.def("basis", &Lielab::domain::rn::basis);
    Lielab_domain_rn.def("project", &Lielab::domain::rn::project);
    Lielab_domain_rn.def("get_dimension", &Lielab::domain::rn::get_dimension);
    Lielab_domain_rn.def("get_vector", &Lielab::domain::rn::get_vector);
    Lielab_domain_rn.def("get_matrix", &Lielab::domain::rn::get_matrix);
    Lielab_domain_rn.def("set_vector", &Lielab::domain::rn::set_vector);
    Lielab_domain_rn.def("__call__", [](const Lielab::domain::rn & self, const size_t index)
        {
            return self(index);
        });
    Lielab_domain_rn.def("__call__", [](const Lielab::domain::rn & self, const size_t index1, const size_t index2)
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
    Lielab_domain_rn.def(py::self *= py::self);
    Lielab_domain_rn.def(py::self / int());
    Lielab_domain_rn.def(py::self / double());
    Lielab_domain_rn.def(py::self /= int());
    Lielab_domain_rn.def(py::self /= double());
    Lielab_domain_rn.def(int() * py::self);
    Lielab_domain_rn.def(double() * py::self);
    Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::gl  & rhs) {return lhs+rhs;});
    Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    Lielab_domain_rn.def("__add__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::GL  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__mul__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::gl  & rhs) {return lhs-rhs;});
    Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    Lielab_domain_rn.def("__sub__", [](const Lielab::domain::rn  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});
    Lielab_domain_rn.def("__repr__", [](const Lielab::domain::rn & self)
        {
            return "<lielab.domain.rn>";
        });
    Lielab_domain_rn.def("__str__", [](const Lielab::domain::rn & self)
        {
            return matstr(self._data);
        });
    
    /*!
    * Bindings for Lielab::domain::se
    */

    auto Lielab_domain_se = py::class_<Lielab::domain::se>(m_domain, "se");
    Lielab_domain_se.def_readwrite("_data", &Lielab::domain::se::_data);
    Lielab_domain_se.def_readwrite("shape", &Lielab::domain::se::shape);
    Lielab_domain_se.def_readonly_static("abelian", &Lielab::domain::se::abelian);
    Lielab_domain_se.def(py::init<>());
    Lielab_domain_se.def(py::init<const size_t>());
    Lielab_domain_se.def(py::init<const Eigen::MatrixXd &>());
    Lielab_domain_se.def("basis", &Lielab::domain::se::basis);
        // .def("project", &Lielab::domain::se::project);
    Lielab_domain_se.def("get_dimension", &Lielab::domain::se::get_dimension);
    Lielab_domain_se.def("get_vector", &Lielab::domain::se::get_vector);
    Lielab_domain_se.def("get_matrix", &Lielab::domain::se::get_matrix);
    Lielab_domain_se.def("set_vector", &Lielab::domain::se::set_vector);
    Lielab_domain_se.def("__call__", [](const Lielab::domain::se & self, const size_t index)
        {
            return self(index);
        });
    Lielab_domain_se.def("__call__", [](const Lielab::domain::se & self, const size_t index1, const size_t index2)
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
    Lielab_domain_se.def(py::self *= py::self);
    Lielab_domain_se.def(py::self / int());
    Lielab_domain_se.def(py::self / double());
    Lielab_domain_se.def(py::self /= int());
    Lielab_domain_se.def(py::self /= double());
    Lielab_domain_se.def(int() * py::self);
    Lielab_domain_se.def(double() * py::self);
    Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::gl  & rhs) {return lhs+rhs;});
    Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    Lielab_domain_se.def("__add__", [](const Lielab::domain::se  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::GL  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__mul__", [](const Lielab::domain::se  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::gl  & rhs) {return lhs-rhs;});
    Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    Lielab_domain_se.def("__sub__", [](const Lielab::domain::se  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});
    Lielab_domain_se.def("__repr__", [](const Lielab::domain::se & self)
        {
            return "<lielab.domain.se>";
        });
    Lielab_domain_se.def("__str__", [](const Lielab::domain::se & self)
        {
            return matstr(self._data);
        });


    /*!
    * Bindings for Lielab::domain::so
    */

    auto Lielab_domain_so = py::class_<Lielab::domain::so>(m_domain, "so");
    Lielab_domain_so.def_readwrite("_data", &Lielab::domain::so::_data);
    Lielab_domain_so.def_readwrite("shape", &Lielab::domain::so::shape);
    Lielab_domain_so.def_readonly_static("abelian", &Lielab::domain::so::abelian);
    Lielab_domain_so.def(py::init<>());
    Lielab_domain_so.def(py::init<const size_t>());
    Lielab_domain_so.def(py::init<const Eigen::MatrixXd &>());
    Lielab_domain_so.def("basis", &Lielab::domain::so::basis);
    Lielab_domain_so.def("project", &Lielab::domain::so::project);
    Lielab_domain_so.def("get_dimension", &Lielab::domain::so::get_dimension);
    Lielab_domain_so.def("get_vector", &Lielab::domain::so::get_vector);
    Lielab_domain_so.def("get_matrix", &Lielab::domain::so::get_matrix);
    Lielab_domain_so.def("set_vector", &Lielab::domain::so::set_vector);
    Lielab_domain_so.def("__call__", [](const Lielab::domain::so & self, const size_t index)
        {
            return self(index);
        });
    Lielab_domain_so.def("__call__", [](const Lielab::domain::so & self, const size_t index1, const size_t index2)
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
    Lielab_domain_so.def(py::self *= py::self);
    Lielab_domain_so.def(py::self / int());
    Lielab_domain_so.def(py::self / double());
    Lielab_domain_so.def(py::self /= int());
    Lielab_domain_so.def(py::self /= double());
    Lielab_domain_so.def(int() * py::self);
    Lielab_domain_so.def(double() * py::self);
    Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::gl  & rhs) {return lhs+rhs;});
    Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    Lielab_domain_so.def("__add__", [](const Lielab::domain::so  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::GL  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__mul__", [](const Lielab::domain::so  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::gl  & rhs) {return lhs-rhs;});
    Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    Lielab_domain_so.def("__sub__", [](const Lielab::domain::so  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});
    Lielab_domain_so.def("__repr__", [](const Lielab::domain::so & self)
        {
            return "<lielab.domain.so>";
        });
    Lielab_domain_so.def("__str__", [](const Lielab::domain::so & self)
        {
            return matstr(self._data);
        });

    
    /*!
    * Bindings for Lielab::domain::sp
    */
    
    auto Lielab_domain_sp = py::class_<Lielab::domain::sp>(m_domain, "sp");
    Lielab_domain_sp.def_readwrite("_data", &Lielab::domain::sp::_data);
    Lielab_domain_sp.def_readwrite("shape", &Lielab::domain::sp::shape);
    Lielab_domain_sp.def_readonly_static("abelian", &Lielab::domain::sp::abelian);
    Lielab_domain_sp.def(py::init<>());
    Lielab_domain_sp.def(py::init<const size_t>());
    Lielab_domain_sp.def(py::init<const Eigen::MatrixXd &>());
    Lielab_domain_sp.def("basis", &Lielab::domain::sp::basis);
    Lielab_domain_sp.def("project", &Lielab::domain::sp::project);
    Lielab_domain_sp.def("get_dimension", &Lielab::domain::sp::get_dimension);
    Lielab_domain_sp.def("get_vector", &Lielab::domain::sp::get_vector);
    Lielab_domain_sp.def("get_matrix", &Lielab::domain::sp::get_matrix);
    Lielab_domain_sp.def("set_vector", &Lielab::domain::sp::set_vector);
    Lielab_domain_sp.def("__call__", [](const Lielab::domain::sp & self, const size_t index)
        {
            return self(index);
        });
    Lielab_domain_sp.def("__call__", [](const Lielab::domain::sp & self, const size_t index1, const size_t index2)
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
    Lielab_domain_sp.def(py::self *= py::self);
    Lielab_domain_sp.def(py::self / int());
    Lielab_domain_sp.def(py::self / double());
    Lielab_domain_sp.def(py::self /= int());
    Lielab_domain_sp.def(py::self /= double());
    Lielab_domain_sp.def(int() * py::self);
    Lielab_domain_sp.def(double() * py::self);
    Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::gl  & rhs) {return lhs+rhs;});
    Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    Lielab_domain_sp.def("__add__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::GL  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__mul__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::gl  & rhs) {return lhs-rhs;});
    Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    Lielab_domain_sp.def("__sub__", [](const Lielab::domain::sp  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});
    Lielab_domain_sp.def("__repr__", [](const Lielab::domain::sp & self)
        {
            return "<lielab.domain.sp>";
        });
    Lielab_domain_sp.def("__str__", [](const Lielab::domain::sp & self)
        {
            return matstr(self._data);
        });
    
    
    /*!
    * Bindings for Lielab::domain::su
    */
    
    auto Lielab_domain_su = py::class_<Lielab::domain::su>(m_domain, "su");
    Lielab_domain_su.def_readwrite("_data", &Lielab::domain::su::_data);
    Lielab_domain_su.def_readwrite("shape", &Lielab::domain::su::shape);
    Lielab_domain_su.def_readonly_static("abelian", &Lielab::domain::su::abelian);
    Lielab_domain_su.def(py::init<>());
    Lielab_domain_su.def(py::init<const size_t>());
    Lielab_domain_su.def(py::init<const Eigen::MatrixXcd &>());
    Lielab_domain_su.def("basis", &Lielab::domain::su::basis);
        // .def("project", &Lielab::domain::su::project); // TODO: Implement project
    Lielab_domain_su.def("get_dimension", &Lielab::domain::su::get_dimension);
    Lielab_domain_su.def("get_vector", &Lielab::domain::su::get_vector);
    Lielab_domain_su.def("get_matrix", &Lielab::domain::su::get_matrix);
    Lielab_domain_su.def("set_vector", &Lielab::domain::su::set_vector);
    Lielab_domain_su.def("__call__", [](const Lielab::domain::su & self, const size_t index)
        {
            return self(index);
        });
    Lielab_domain_su.def("__call__", [](const Lielab::domain::su & self, const size_t index1, const size_t index2)
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
    Lielab_domain_su.def(py::self *= py::self);
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
    Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::cn  & rhs) {return lhs+rhs;});
    Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::gl  & rhs) {return lhs+rhs;});
    Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::glc & rhs) {return lhs+rhs;});
    Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::rn  & rhs) {return lhs+rhs;});
    Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::se  & rhs) {return lhs+rhs;});
    Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::so  & rhs) {return lhs+rhs;});
    Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::sp  & rhs) {return lhs+rhs;});
    Lielab_domain_su.def("__add__", [](const Lielab::domain::su  & lhs, const Lielab::domain::su  & rhs) {return lhs+rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::CN  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::GL  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::GLC & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::RN  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::SE  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::SO  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::SP  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__mul__", [](const Lielab::domain::su  & lhs, const Lielab::domain::SU  & rhs) {return lhs*rhs;});
    Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::cn  & rhs) {return lhs-rhs;});
    Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::gl  & rhs) {return lhs-rhs;});
    Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::glc & rhs) {return lhs-rhs;});
    Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::rn  & rhs) {return lhs-rhs;});
    Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::se  & rhs) {return lhs-rhs;});
    Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::so  & rhs) {return lhs-rhs;});
    Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::sp  & rhs) {return lhs-rhs;});
    Lielab_domain_su.def("__sub__", [](const Lielab::domain::su  & lhs, const Lielab::domain::su  & rhs) {return lhs-rhs;});
    Lielab_domain_su.def("__repr__", [](const Lielab::domain::su & self)
        {
            return "<lielab.domain.su>";
        });
    Lielab_domain_su.def("__str__", [](const Lielab::domain::su & self)
        {
            return matstr(self._data);
        });
    
    /*!
     * Lie Groups
     */

    auto Lielab_domain_CN = py::class_<Lielab::domain::CN>(m_domain, "CN");
    Lielab_domain_CN.def_readwrite("_data", &Lielab::domain::CN::_data);
    Lielab_domain_CN.def_readwrite("shape", &Lielab::domain::CN::shape);
    Lielab_domain_CN.def_readonly_static("abelian", &Lielab::domain::CN::abelian);
    Lielab_domain_CN.def(py::init<>());
    Lielab_domain_CN.def(py::init<int>());
    Lielab_domain_CN.def(py::init<Eigen::MatrixXcd>());
        // .def("project", &Lielab::domain::CN::project); // TODO:
    Lielab_domain_CN.def("get_dimension", &Lielab::domain::CN::get_dimension);
    Lielab_domain_CN.def("get_matrix", &Lielab::domain::CN::get_matrix);
    Lielab_domain_CN.def("inverse", &Lielab::domain::CN::inverse);
    Lielab_domain_CN.def("serialize", &Lielab::domain::CN::serialize);
    Lielab_domain_CN.def("unserialize", &Lielab::domain::CN::unserialize);
    Lielab_domain_CN.def("__call__", [](const Lielab::domain::CN & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_CN.def(py::self * py::self);
    Lielab_domain_CN.def(py::self *= py::self);
    Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_CN.def("__mul__", [](const Lielab::domain::CN  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_CN.def("__repr__", [](const Lielab::domain::CN & self)
        {
            return "<lielab.domain.CN>";
        });
    Lielab_domain_CN.def("__str__", [](const Lielab::domain::CN & self)
        {
            return matstr(self._data);
        });
    
    auto Lielab_domain_GL = py::class_<Lielab::domain::GL>(m_domain, "GL");
    Lielab_domain_GL.def_readwrite("_data", &Lielab::domain::GL::_data);
    Lielab_domain_GL.def_readonly_static("abelian", &Lielab::domain::GL::abelian);
    Lielab_domain_GL.def(py::init<>());
    Lielab_domain_GL.def(py::init<int>());
    Lielab_domain_GL.def(py::init<Eigen::MatrixXd>());
    Lielab_domain_GL.def("project", &Lielab::domain::GL::project);
    Lielab_domain_GL.def("get_dimension", &Lielab::domain::GL::get_dimension);
    Lielab_domain_GL.def("get_matrix", &Lielab::domain::GL::get_matrix);
    Lielab_domain_GL.def("inverse", &Lielab::domain::GL::inverse);
    Lielab_domain_GL.def("serialize", &Lielab::domain::GL::serialize);
    Lielab_domain_GL.def("unserialize", &Lielab::domain::GL::unserialize);
    Lielab_domain_GL.def("__call__", [](const Lielab::domain::GL & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_GL.def(py::self * py::self);
    Lielab_domain_GL.def(py::self *= py::self);
    Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GL  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GL  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GL  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GL  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GL  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GL  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GL  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_GL.def("__mul__", [](const Lielab::domain::GL  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_GL.def("__repr__", [](const Lielab::domain::GL & self)
        {
            return "<lielab.domain.GL>";
        });
    Lielab_domain_GL.def("__str__", [](const Lielab::domain::GL & self)
        {
            return matstr(self._data);
        });
    
    auto Lielab_domain_GLC = py::class_<Lielab::domain::GLC>(m_domain, "GLC");
    Lielab_domain_GLC.def_readwrite("_data", &Lielab::domain::GLC::_data);
    Lielab_domain_GLC.def_readonly_static("abelian", &Lielab::domain::GLC::abelian);
    Lielab_domain_GLC.def(py::init<>());
    Lielab_domain_GLC.def(py::init<int>());
    Lielab_domain_GLC.def(py::init<Eigen::MatrixXd>());
    Lielab_domain_GLC.def("project", &Lielab::domain::GLC::project);
    Lielab_domain_GLC.def("get_dimension", &Lielab::domain::GLC::get_dimension);
    Lielab_domain_GLC.def("get_matrix", &Lielab::domain::GLC::get_matrix);
    Lielab_domain_GLC.def("inverse", &Lielab::domain::GLC::inverse);
    Lielab_domain_GLC.def("serialize", &Lielab::domain::GLC::serialize);
    Lielab_domain_GLC.def("unserialize", &Lielab::domain::GLC::unserialize);
    Lielab_domain_GLC.def("__call__", [](const Lielab::domain::GLC & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_GLC.def(py::self * py::self);
    Lielab_domain_GLC.def(py::self *= py::self);
    Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_GLC.def("__mul__", [](const Lielab::domain::GLC & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_GLC.def("__repr__", [](const Lielab::domain::GLC & self)
        {
            return "<lielab.domain.GLC>";
        });
    Lielab_domain_GLC.def("__str__", [](const Lielab::domain::GLC & self)
        {
            return matstr(self._data);
        });

    auto Lielab_domain_RN = py::class_<Lielab::domain::RN>(m_domain, "RN");
    Lielab_domain_RN.def_readwrite("_data", &Lielab::domain::RN::_data);
    Lielab_domain_RN.def_readwrite("shape", &Lielab::domain::RN::shape);
    Lielab_domain_RN.def_readonly_static("abelian", &Lielab::domain::RN::abelian);
    Lielab_domain_RN.def(py::init<>());
    Lielab_domain_RN.def(py::init<int>());
    Lielab_domain_RN.def(py::init<Eigen::MatrixXd>());
    Lielab_domain_RN.def("project", &Lielab::domain::RN::project);
    Lielab_domain_RN.def("get_dimension", &Lielab::domain::RN::get_dimension);
    Lielab_domain_RN.def("get_matrix", &Lielab::domain::RN::get_matrix);
    Lielab_domain_RN.def("inverse", &Lielab::domain::RN::inverse);
    Lielab_domain_RN.def("serialize", &Lielab::domain::RN::serialize);
    Lielab_domain_RN.def("unserialize", &Lielab::domain::RN::unserialize);
    Lielab_domain_RN.def("__call__", [](const Lielab::domain::RN & self, const size_t index)
        {
            return self(index);
        });
    Lielab_domain_RN.def("__call__", [](const Lielab::domain::RN & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_RN.def(py::self * py::self);
    Lielab_domain_RN.def(py::self *= py::self);
    Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_RN.def("__mul__", [](const Lielab::domain::RN  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_RN.def("__repr__", [](const Lielab::domain::RN & self)
        {
            return "<lielab.domain.RN>";
        });
    Lielab_domain_RN.def("__str__", [](const Lielab::domain::RN & self)
        {
            return matstr(self._data);
        });

    auto Lielab_domain_SE = py::class_<Lielab::domain::SE>(m_domain, "SE");
    Lielab_domain_SE.def_readwrite("_data", &Lielab::domain::SE::_data);
    Lielab_domain_SE.def_readwrite("shape", &Lielab::domain::SE::shape);
    Lielab_domain_SE.def_readonly_static("abelian", &Lielab::domain::SE::abelian);
    Lielab_domain_SE.def(py::init<>());
    Lielab_domain_SE.def(py::init<int>());
    Lielab_domain_SE.def(py::init<Eigen::MatrixXd>());
        // .def("project", &Lielab::domain::SE::project);
    Lielab_domain_SE.def("get_dimension", &Lielab::domain::SE::get_dimension);
    Lielab_domain_SE.def("get_matrix", &Lielab::domain::SE::get_matrix);
    Lielab_domain_SE.def("inverse", &Lielab::domain::SE::inverse);
    Lielab_domain_SE.def("serialize", &Lielab::domain::SE::serialize);
    Lielab_domain_SE.def("unserialize", &Lielab::domain::SE::unserialize);
    Lielab_domain_SE.def("__call__", [](const Lielab::domain::SE & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SE.def(py::self * py::self);
    Lielab_domain_SE.def(py::self *= py::self);
    Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_SE.def("__mul__", [](const Lielab::domain::SE  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_SE.def("__repr__", [](const Lielab::domain::SE & self)
        {
            return "<lielab.domain.SE>";
        });
    Lielab_domain_SE.def("__str__", [](const Lielab::domain::SE & self)
        {
            return matstr(self._data);
        });

    auto Lielab_domain_SO = py::class_<Lielab::domain::SO>(m_domain, "SO");
    Lielab_domain_SO.def_readwrite("_data", &Lielab::domain::SO::_data);
    Lielab_domain_SO.def_readwrite("shape", &Lielab::domain::SO::shape);
    Lielab_domain_SO.def_readonly_static("abelian", &Lielab::domain::SO::abelian);
    Lielab_domain_SO.def(py::init<>());
    Lielab_domain_SO.def(py::init<int>());
    Lielab_domain_SO.def(py::init<Eigen::MatrixXd>());
    Lielab_domain_SO.def("project", &Lielab::domain::SO::project);
    Lielab_domain_SO.def("get_dimension", &Lielab::domain::SO::get_dimension);
    Lielab_domain_SO.def("get_matrix", &Lielab::domain::SO::get_matrix);
    Lielab_domain_SO.def("inverse", &Lielab::domain::SO::inverse);
    Lielab_domain_SO.def("serialize", &Lielab::domain::SO::serialize);
    Lielab_domain_SO.def("unserialize", &Lielab::domain::SO::unserialize);
    Lielab_domain_SO.def("__call__", [](const Lielab::domain::SO & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SO.def(py::self * py::self);
    Lielab_domain_SO.def(py::self *= py::self);
    Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_SO.def("__mul__", [](const Lielab::domain::SO  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_SO.def("__repr__", [](const Lielab::domain::SO & self)
        {
            return "<lielab.domain.SO>";
        });
    Lielab_domain_SO.def("__str__", [](const Lielab::domain::SO & self)
        {
            return matstr(self._data);
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
    

    auto Lielab_domain_SP = py::class_<Lielab::domain::SP>(m_domain, "SP");
    Lielab_domain_SP.def_readwrite("_data", &Lielab::domain::SP::_data);
    Lielab_domain_SP.def_readwrite("shape", &Lielab::domain::SP::shape);
    Lielab_domain_SP.def_readonly_static("abelian", &Lielab::domain::SP::abelian);
    Lielab_domain_SP.def(py::init<>());
    Lielab_domain_SP.def(py::init<int>());
    Lielab_domain_SP.def(py::init<Eigen::MatrixXd>());
        // .def("project", &Lielab::domain::SP::project); // TODO:
    Lielab_domain_SP.def("get_dimension", &Lielab::domain::SP::get_dimension);
    Lielab_domain_SP.def("get_matrix", &Lielab::domain::SP::get_matrix);
    Lielab_domain_SP.def("inverse", &Lielab::domain::SP::inverse);
    Lielab_domain_SP.def("serialize", &Lielab::domain::SP::serialize);
    Lielab_domain_SP.def("unserialize", &Lielab::domain::SP::unserialize);
    Lielab_domain_SP.def("__call__", [](const Lielab::domain::SP & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SP.def(py::self * py::self);
    Lielab_domain_SP.def(py::self *= py::self);
    Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_SP.def("__mul__", [](const Lielab::domain::SP  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_SP.def("__repr__", [](const Lielab::domain::SP & self)
        {
            return "<lielab.domain.SP>";
        });
    Lielab_domain_SP.def("__str__", [](const Lielab::domain::SP & self)
        {
            return matstr(self._data);
        });

    auto Lielab_domain_SU = py::class_<Lielab::domain::SU>(m_domain, "SU");
    Lielab_domain_SU.def_readwrite("_data", &Lielab::domain::SU::_data);
    Lielab_domain_SU.def_readwrite("shape", &Lielab::domain::SU::shape);
    Lielab_domain_SU.def_readonly_static("abelian", &Lielab::domain::SU::abelian);
    Lielab_domain_SU.def(py::init<>());
    Lielab_domain_SU.def(py::init<int>());
    Lielab_domain_SU.def(py::init<Eigen::MatrixXcd>());
        // .def("project", &Lielab::domain::SU::project); // TODO:
    Lielab_domain_SU.def("get_dimension", &Lielab::domain::SU::get_dimension);
    Lielab_domain_SU.def("get_matrix", &Lielab::domain::SU::get_matrix);
    Lielab_domain_SU.def("inverse", &Lielab::domain::SU::inverse);
    Lielab_domain_SU.def("serialize", &Lielab::domain::SU::serialize);
    Lielab_domain_SU.def("unserialize", &Lielab::domain::SU::unserialize);
    Lielab_domain_SU.def("__call__", [](const Lielab::domain::SU & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_SU.def(py::self * py::self);
    Lielab_domain_SU.def(py::self *= py::self);
    Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::cn  & rhs) {return lhs*rhs;});
    Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::gl  & rhs) {return lhs*rhs;});
    Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::glc & rhs) {return lhs*rhs;});
    Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::rn  & rhs) {return lhs*rhs;});
    Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::se  & rhs) {return lhs*rhs;});
    Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::so  & rhs) {return lhs*rhs;});
    Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::sp  & rhs) {return lhs*rhs;});
    Lielab_domain_SU.def("__mul__", [](const Lielab::domain::SU  & lhs, const Lielab::domain::su  & rhs) {return lhs*rhs;});
    Lielab_domain_SU.def("__repr__", [](const Lielab::domain::SU & self)
        {
            return "<lielab.domain.SU>";
        });
    Lielab_domain_SU.def("__str__", [](const Lielab::domain::SU & self)
        {
            return matstr(self._data);
        });
    Lielab_domain_SU.def_static("from_quaternion", &Lielab::domain::SU::from_quaternion<double>);

    auto Lielab_domain_CompositeAlgebra = py::class_<Lielab::domain::CompositeAlgebra>(m_domain, "CompositeAlgebra");
    Lielab_domain_CompositeAlgebra.def(py::init());
    Lielab_domain_CompositeAlgebra.def(py::init<std::vector<Lielab::domain::CompositeAlgebra::TYPES>>());
    Lielab_domain_CompositeAlgebra.def_readwrite("space", &Lielab::domain::CompositeAlgebra::space);
    Lielab_domain_CompositeAlgebra.def("get_dimension", &Lielab::domain::CompositeAlgebra::get_dimension);
    Lielab_domain_CompositeAlgebra.def("get_shape", &Lielab::domain::CompositeAlgebra::get_shape);
    Lielab_domain_CompositeAlgebra.def("get_vector", &Lielab::domain::CompositeAlgebra::get_vector);
    Lielab_domain_CompositeAlgebra.def("set_vector", &Lielab::domain::CompositeAlgebra::set_vector);
    Lielab_domain_CompositeAlgebra.def("get_matrix", &Lielab::domain::CompositeAlgebra::get_matrix);
    Lielab_domain_CompositeAlgebra.def(py::self + py::self);
    Lielab_domain_CompositeAlgebra.def(py::self += py::self);
    Lielab_domain_CompositeAlgebra.def(py::self - py::self);
    Lielab_domain_CompositeAlgebra.def(py::self -= py::self);
    Lielab_domain_CompositeAlgebra.def(-py::self);
    Lielab_domain_CompositeAlgebra.def(py::self * int());
    Lielab_domain_CompositeAlgebra.def(py::self * double());
    Lielab_domain_CompositeAlgebra.def(py::self * py::self);
    Lielab_domain_CompositeAlgebra.def(py::self *= int());
    Lielab_domain_CompositeAlgebra.def(py::self *= double());
    Lielab_domain_CompositeAlgebra.def(py::self *= py::self);
    Lielab_domain_CompositeAlgebra.def(py::self / int());
    Lielab_domain_CompositeAlgebra.def(py::self / double());
    Lielab_domain_CompositeAlgebra.def(py::self /= int());
    Lielab_domain_CompositeAlgebra.def(py::self /= double());
    Lielab_domain_CompositeAlgebra.def(int() * py::self);
    Lielab_domain_CompositeAlgebra.def(double() * py::self);
    Lielab_domain_CompositeAlgebra.def("__repr__",
        [](const Lielab::domain::CompositeAlgebra & self)
        {
            return "<lielab.domain.CompositeAlgebra>";
        });
    Lielab_domain_CompositeAlgebra.def("__str__",
        [](const Lielab::domain::CompositeAlgebra & self)
        {
            return self.to_string();
        });
    

    auto Lielab_domain_CompositeManifold = py::class_<Lielab::domain::CompositeManifold>(m_domain, "CompositeManifold");
    Lielab_domain_CompositeManifold.def(py::init());
    Lielab_domain_CompositeManifold.def(py::init<std::vector<Lielab::domain::CompositeManifold::TYPES>>());
    Lielab_domain_CompositeManifold.def_readwrite("space", &Lielab::domain::CompositeManifold::space);
    Lielab_domain_CompositeManifold.def("get_shape", &Lielab::domain::CompositeManifold::get_shape);
    Lielab_domain_CompositeManifold.def("serialize", &Lielab::domain::CompositeManifold::serialize);
    Lielab_domain_CompositeManifold.def("unserialize", &Lielab::domain::CompositeManifold::unserialize);
    Lielab_domain_CompositeManifold.def("get_matrix", &Lielab::domain::CompositeManifold::get_matrix);
    Lielab_domain_CompositeManifold.def(py::self * py::self);
    Lielab_domain_CompositeManifold.def(py::self *= py::self);
    Lielab_domain_CompositeManifold.def("inverse", &Lielab::domain::CompositeManifold::inverse);
    Lielab_domain_CompositeManifold.def("__repr__",
        [](const Lielab::domain::CompositeManifold & self)
        {
            return "<lielab.domain.CompositeManifold>";
        });
    Lielab_domain_CompositeManifold.def("__str__",
        [](const Lielab::domain::CompositeManifold & self)
        {
            return self.to_string();
        });


    /*!
    * Begin content for the "functions" submodule.
    */

    py::module m_functions = m.def_submodule("functions", "The functions submodule.");
    // m_functions.def("copair", &Lielab::functions::copair<Lielab::domain::gl>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<Lielab::domain::rn>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<Lielab::domain::so>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<Lielab::domain::sp>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<Lielab::domain::su>, "The copair function.");
    m_functions.def("factorial", &Lielab::functions::factorial, "The factorial function");
    m_functions.def("Ad_numerical", [](const Lielab::domain::cn & a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"), "The numerical Ad function.");
    m_functions.def("Ad_numerical", [](const Lielab::domain::gl & a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"), "The numerical Ad function.");
    m_functions.def("Ad_numerical", [](const Lielab::domain::glc & a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"), "The numerical Ad function.");
    m_functions.def("Ad_numerical", [](const Lielab::domain::rn & a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"), "The numerical Ad function.");
    m_functions.def("Ad_numerical", [](const Lielab::domain::se & a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"), "The numerical Ad function.");
    m_functions.def("Ad_numerical", [](const Lielab::domain::so & a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"), "The numerical Ad function.");
    m_functions.def("Ad_numerical", [](const Lielab::domain::sp & a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"), "The numerical Ad function.");
    m_functions.def("Ad_numerical", [](const Lielab::domain::su & a){return Lielab::functions::Ad_numerical(a);}, py::arg("a"), "The numerical Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::cn & a){return Lielab::functions::Ad(a);}, py::arg("a"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::gl & a){return Lielab::functions::Ad(a);}, py::arg("a"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::glc & a){return Lielab::functions::Ad(a);}, py::arg("a"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::rn & a){return Lielab::functions::Ad(a);}, py::arg("a"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::se & a){return Lielab::functions::Ad(a);}, py::arg("a"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::so & a){return Lielab::functions::Ad(a);}, py::arg("a"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::sp & a){return Lielab::functions::Ad(a);}, py::arg("a"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::su & a){return Lielab::functions::Ad(a);}, py::arg("a"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::cn & a, Lielab::domain::cn & b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::gl & a, Lielab::domain::gl & b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::glc & a, Lielab::domain::glc & b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::rn & a, Lielab::domain::rn & b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::se & a, Lielab::domain::se & b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::so & a, Lielab::domain::so & b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::sp & a, Lielab::domain::sp & b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::su & a, Lielab::domain::su & b){return Lielab::functions::Ad(a, b);}, py::arg("a"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::CN & A, Lielab::domain::cn & b){return Lielab::functions::Ad(A, b);}, py::arg("A"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::GL & A, Lielab::domain::gl & b){return Lielab::functions::Ad(A, b);}, py::arg("A"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::GLC & A, Lielab::domain::glc & b){return Lielab::functions::Ad(A, b);}, py::arg("A"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::RN & A, Lielab::domain::rn & b){return Lielab::functions::Ad(A, b);}, py::arg("A"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::SE & A, Lielab::domain::se & b){return Lielab::functions::Ad(A, b);}, py::arg("A"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::SO & A, Lielab::domain::so & b){return Lielab::functions::Ad(A, b);}, py::arg("A"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::SP & A, Lielab::domain::sp & b){return Lielab::functions::Ad(A, b);}, py::arg("A"), py::arg("b"), "The Ad function.");
    m_functions.def("Ad", [](const Lielab::domain::SU & A, Lielab::domain::su & b){return Lielab::functions::Ad(A, b);}, py::arg("A"), py::arg("b"), "The Ad function.");
    // m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::cn>, "The cayley1 function.");
    // m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::gl>, "The cayley1 function.");
    // m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::glc>, "The cayley1 function.");
    m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::rn>, "The cayley1 function.");
    // m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::se>, "The cayley1 function.");
    m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::so>, "The cayley1 function.");
    m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::sp>, "The cayley1 function.");
    m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::su>, "The cayley1 function.");
    // m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::cn>, "The cayley2 function.");
    // m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::gl>, "The cayley2 function.");
    // m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::glc>, "The cayley2 function.");
    m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::rn>, "The cayley2 function.");
    // m_functions.def("cayley2", &Lielab::functions::cayley1<Lielab::domain::se>, "The cayley2 function.");
    m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::so>, "The cayley2 function.");
    m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::sp>, "The cayley2 function.");
    // m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::su>, "The cayley2 function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::cn>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::gl>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::glc>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::rn>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::se>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::so>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::sp>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::su>, "The commutator function.");
    // m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::cn>, "The Killing function.");
    // m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::gl>, "The Killing function.");
    // m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::glc>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::rn>, "The Killing function.");
    // m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::se>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::so>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::sp>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::su>, "The Killing function."); // TODO: Might be wrong
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::cn>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::gl>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::glc>, "The Killingform function.");
    // m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::se>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::rn>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::so>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::sp>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::su>, "The Killingform function."); // TODO: Might be wrong
    m_functions.def("ad_numerical", [](const Lielab::domain::cn & a, const int p){return Lielab::functions::ad_numerical(a, p);}, py::arg("a"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::cn & a, const Lielab::domain::cn & b, const int p){return Lielab::functions::ad_numerical(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::gl & a, const int p){return Lielab::functions::ad_numerical(a, p);}, py::arg("a"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::gl & a, const Lielab::domain::gl & b, const int p){return Lielab::functions::ad_numerical(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::glc & a, const int p){return Lielab::functions::ad_numerical(a, p);}, py::arg("a"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::glc & a, const Lielab::domain::glc & b, const int p){return Lielab::functions::ad_numerical(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::rn & a, const int p){return Lielab::functions::ad_numerical(a, p);}, py::arg("a"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::rn & a, const Lielab::domain::rn & b, const int p){return Lielab::functions::ad_numerical(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::se & a, const int p){return Lielab::functions::ad_numerical(a, p);}, py::arg("a"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::se & a, const Lielab::domain::se & b, const int p){return Lielab::functions::ad_numerical(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::so & a, const int p){return Lielab::functions::ad_numerical(a, p);}, py::arg("a"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::so & a, const Lielab::domain::so & b, const int p){return Lielab::functions::ad_numerical(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::sp & a, const int p){return Lielab::functions::ad_numerical(a, p);}, py::arg("a"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::sp & a, const Lielab::domain::sp & b, const int p){return Lielab::functions::ad_numerical(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::su & a, const int p){return Lielab::functions::ad_numerical(a, p);}, py::arg("a"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad_numerical", [](const Lielab::domain::su & a, const Lielab::domain::su & b, const int p){return Lielab::functions::ad_numerical(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The numerical ad function.");
    m_functions.def("ad", [](const Lielab::domain::cn & a, const int p){return Lielab::functions::ad(a, p);}, py::arg("a"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::cn & a, const Lielab::domain::cn & b, const int p){return Lielab::functions::ad(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::gl & a, const int p){return Lielab::functions::ad(a, p);}, py::arg("a"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::gl & a, const Lielab::domain::gl & b, const int p){return Lielab::functions::ad(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::glc & a, const int p){return Lielab::functions::ad(a, p);}, py::arg("a"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::glc & a, const Lielab::domain::glc & b, const int p){return Lielab::functions::ad(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::rn & a, const int p){return Lielab::functions::ad(a, p);}, py::arg("a"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::rn & a, const Lielab::domain::rn & b, const int p){return Lielab::functions::ad(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::se & a, const int p){return Lielab::functions::ad(a, p);}, py::arg("a"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::se & a, const Lielab::domain::se & b, const int p){return Lielab::functions::ad(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::so & a, const int p){return Lielab::functions::ad(a, p);}, py::arg("a"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::so & a, const Lielab::domain::so & b, const int p){return Lielab::functions::ad(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::sp & a, const int p){return Lielab::functions::ad(a, p);}, py::arg("a"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::sp & a, const Lielab::domain::sp & b, const int p){return Lielab::functions::ad(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::su & a, const int p){return Lielab::functions::ad(a, p);}, py::arg("a"), py::arg("p") = 1, "The ad function.");
    m_functions.def("ad", [](const Lielab::domain::su & a, const Lielab::domain::su & b, const int p){return Lielab::functions::ad(a, b, p);}, py::arg("a"), py::arg("b"), py::arg("p") = 1, "The ad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::cn & a){return Lielab::functions::coad_numerical(a);}, py::arg("a"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::cn & a, const Lielab::domain::cn & b){return Lielab::functions::coad_numerical(a, b);}, py::arg("a"), py::arg("b"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::gl & a){return Lielab::functions::coad_numerical(a);}, py::arg("a"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::gl & a, const Lielab::domain::gl & b){return Lielab::functions::coad_numerical(a, b);}, py::arg("a"), py::arg("b"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::glc & a){return Lielab::functions::coad_numerical(a);}, py::arg("a"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::glc & a, const Lielab::domain::glc & b){return Lielab::functions::coad_numerical(a, b);}, py::arg("a"), py::arg("b"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::rn & a){return Lielab::functions::coad_numerical(a);}, py::arg("a"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::rn & a, const Lielab::domain::rn & b){return Lielab::functions::coad_numerical(a, b);}, py::arg("a"), py::arg("b"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::se & a){return Lielab::functions::coad_numerical(a);}, py::arg("a"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::se & a, const Lielab::domain::se & b){return Lielab::functions::coad_numerical(a, b);}, py::arg("a"), py::arg("b"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::so & a){return Lielab::functions::coad_numerical(a);}, py::arg("a"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::so & a, const Lielab::domain::so & b){return Lielab::functions::coad_numerical(a, b);}, py::arg("a"), py::arg("b"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::sp & a){return Lielab::functions::coad_numerical(a);}, py::arg("a"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::sp & a, const Lielab::domain::sp & b){return Lielab::functions::coad_numerical(a, b);}, py::arg("a"), py::arg("b"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::su & a){return Lielab::functions::coad_numerical(a);}, py::arg("a"), "The numerical coad function.");
    m_functions.def("coad_numerical", [](const Lielab::domain::su & a, const Lielab::domain::su & b){return Lielab::functions::coad_numerical(a, b);}, py::arg("a"), py::arg("b"), "The numerical coad function.");
    m_functions.def("coad", [](const Lielab::domain::cn & a){return Lielab::functions::coad(a);}, py::arg("a"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::cn & a, const Lielab::domain::cn & b){return Lielab::functions::coad(a, b);}, py::arg("a"), py::arg("b"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::gl & a){return Lielab::functions::coad(a);}, py::arg("a"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::gl & a, const Lielab::domain::gl & b){return Lielab::functions::coad(a, b);}, py::arg("a"), py::arg("b"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::glc & a){return Lielab::functions::coad(a);}, py::arg("a"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::glc & a, const Lielab::domain::glc & b){return Lielab::functions::coad(a, b);}, py::arg("a"), py::arg("b"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::rn & a){return Lielab::functions::coad(a);}, py::arg("a"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::rn & a, const Lielab::domain::rn & b){return Lielab::functions::coad(a, b);}, py::arg("a"), py::arg("b"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::se & a){return Lielab::functions::coad(a);}, py::arg("a"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::se & a, const Lielab::domain::se & b){return Lielab::functions::coad(a, b);}, py::arg("a"), py::arg("b"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::so & a){return Lielab::functions::coad(a);}, py::arg("a"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::so & a, const Lielab::domain::so & b){return Lielab::functions::coad(a, b);}, py::arg("a"), py::arg("b"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::sp & a){return Lielab::functions::coad(a);}, py::arg("a"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::sp & a, const Lielab::domain::sp & b){return Lielab::functions::coad(a, b);}, py::arg("a"), py::arg("b"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::su & a){return Lielab::functions::coad(a);}, py::arg("a"), "The coad function.");
    m_functions.def("coad", [](const Lielab::domain::su & a, const Lielab::domain::su & b){return Lielab::functions::coad(a, b);}, py::arg("a"), py::arg("b"), "The coad function.");
    m_functions.def("exp_numerical", &Lielab::functions::exp_numerical<Lielab::domain::cn>, "The numerical exponential function.");
    m_functions.def("exp_numerical", &Lielab::functions::exp_numerical<Lielab::domain::gl>, "The numerical exponential function.");
    m_functions.def("exp_numerical", &Lielab::functions::exp_numerical<Lielab::domain::glc>, "The numerical exponential function.");
    m_functions.def("exp_numerical", &Lielab::functions::exp_numerical<Lielab::domain::rn>, "The numerical exponential function.");
    m_functions.def("exp_numerical", &Lielab::functions::exp_numerical<Lielab::domain::se>, "The numerical exponential function.");
    m_functions.def("exp_numerical", &Lielab::functions::exp_numerical<Lielab::domain::so>, "The numerical exponential function.");
    m_functions.def("exp_numerical", &Lielab::functions::exp_numerical<Lielab::domain::sp>, "The numerical exponential function.");
    m_functions.def("exp_numerical", &Lielab::functions::exp_numerical<Lielab::domain::su>, "The numerical exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::cn>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::gl>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::glc>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::rn>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::se>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::so>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::sp>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::su>, "The exponential function.");
    m_functions.def("log_numerical", &Lielab::functions::log_numerical<Lielab::domain::CN>, "The numerical logarithm function.");
    m_functions.def("log_numerical", &Lielab::functions::log_numerical<Lielab::domain::GL>, "The numerical logarithm function.");
    m_functions.def("log_numerical", &Lielab::functions::log_numerical<Lielab::domain::GLC>, "The numerical logarithm function.");
    m_functions.def("log_numerical", &Lielab::functions::log_numerical<Lielab::domain::RN>, "The numerical logarithm function.");
    m_functions.def("log_numerical", &Lielab::functions::log_numerical<Lielab::domain::SE>, "The numerical logarithm function.");
    m_functions.def("log_numerical", &Lielab::functions::log_numerical<Lielab::domain::SO>, "The numerical logarithm function.");
    m_functions.def("log_numerical", &Lielab::functions::log_numerical<Lielab::domain::SP>, "The numerical logarithm function.");
    m_functions.def("log_numerical", &Lielab::functions::log_numerical<Lielab::domain::SU>, "The numerical logarithm function.");
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::CN>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::GL>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::GLC>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::RN>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::SE>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::SO>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::SP>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::SU>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("bernoulli", &Lielab::functions::bernoulli, "The bernoulli function.");
    // m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::cn>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::gl>, "The dcayley1inv function.");
    // m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::glc>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::rn>, "The dcayley1inv function.");
    // m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::se>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::so>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::sp>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::su>, "The dcayley1inv function.");
    m_functions.def("dexp_numerical", [](const Lielab::domain::cn & a, const size_t order){return Lielab::functions::dexp_numerical(a, order);}, "The numerical dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::cn & a, const Lielab::domain::cn & b, const size_t order){return Lielab::functions::dexp_numerical(a, b, order);}, "The numerical dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::gl & a, const size_t order){return Lielab::functions::dexp_numerical(a, order);}, "The numerical dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::gl & a, const Lielab::domain::gl & b, const size_t order){return Lielab::functions::dexp_numerical(a, b, order);}, "The numerical dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::glc & a, const size_t order){return Lielab::functions::dexp_numerical(a, order);}, "The numerical dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::glc & a, const Lielab::domain::glc & b, const size_t order){return Lielab::functions::dexp_numerical(a, b, order);}, "The numerical dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::rn & a, const size_t order){return Lielab::functions::dexp_numerical(a, order);}, "The numerical dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::rn & a, const Lielab::domain::rn & b, const size_t order){return Lielab::functions::dexp_numerical(a, b, order);}, "The numerical dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::se & a, const size_t order){return Lielab::functions::dexp_numerical(a, order);}, "The numerical dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::se & a, const Lielab::domain::se & b, const size_t order){return Lielab::functions::dexp_numerical(a, b, order);}, "The numerical dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::so & a, const size_t order){return Lielab::functions::dexp_numerical(a, order);}, "The numerical dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::so & a, const Lielab::domain::so & b, const size_t order){return Lielab::functions::dexp_numerical(a, b, order);}, "The numerical dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::sp & a, const size_t order){return Lielab::functions::dexp_numerical(a, order);}, "The numerical dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::sp & a, const Lielab::domain::sp & b, const size_t order){return Lielab::functions::dexp_numerical(a, b, order);}, "The numerical dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::su & a, const size_t order){return Lielab::functions::dexp_numerical(a, order);}, "The numerical dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp_numerical", [](const Lielab::domain::su & a, const Lielab::domain::su & b, const size_t order){return Lielab::functions::dexp_numerical(a, b, order);}, "The numerical dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::cn & a, const size_t order){return Lielab::functions::dexp(a, order);}, "The dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::cn & a, const Lielab::domain::cn & b, const size_t order){return Lielab::functions::dexp(a, b, order);}, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::gl & a, const size_t order){return Lielab::functions::dexp(a, order);}, "The dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::gl & a, const Lielab::domain::gl & b, const size_t order){return Lielab::functions::dexp(a, b, order);}, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::glc & a, const size_t order){return Lielab::functions::dexp(a, order);}, "The dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::glc & a, const Lielab::domain::glc & b, const size_t order){return Lielab::functions::dexp(a, b, order);}, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::rn & a, const size_t order){return Lielab::functions::dexp(a, order);}, "The dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::rn & a, const Lielab::domain::rn & b, const size_t order){return Lielab::functions::dexp(a, b, order);}, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::se & a, const size_t order){return Lielab::functions::dexp(a, order);}, "The dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::se & a, const Lielab::domain::se & b, const size_t order){return Lielab::functions::dexp(a, b, order);}, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::so & a, const size_t order){return Lielab::functions::dexp(a, order);}, "The dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::so & a, const Lielab::domain::so & b, const size_t order){return Lielab::functions::dexp(a, b, order);}, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::sp & a, const size_t order){return Lielab::functions::dexp(a, order);}, "The dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::sp & a, const Lielab::domain::sp & b, const size_t order){return Lielab::functions::dexp(a, b, order);}, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::su & a, const size_t order){return Lielab::functions::dexp(a, order);}, "The dexp function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexp", [](const Lielab::domain::su & a, const Lielab::domain::su & b, const size_t order){return Lielab::functions::dexp(a, b, order);}, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::cn & a, const size_t order){return Lielab::functions::dexpinv_numerical(a, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::cn & a, const Lielab::domain::cn & b, const size_t order){return Lielab::functions::dexpinv_numerical(a, b, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::gl & a, const size_t order){return Lielab::functions::dexpinv_numerical(a, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::gl & a, const Lielab::domain::gl & b, const size_t order){return Lielab::functions::dexpinv_numerical(a, b, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::glc & a, const size_t order){return Lielab::functions::dexpinv_numerical(a, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::glc & a, const Lielab::domain::glc & b, const size_t order){return Lielab::functions::dexpinv_numerical(a, b, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::rn & a, const size_t order){return Lielab::functions::dexpinv_numerical(a, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::rn & a, const Lielab::domain::rn & b, const size_t order){return Lielab::functions::dexpinv_numerical(a, b, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::se & a, const size_t order){return Lielab::functions::dexpinv_numerical(a, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::se & a, const Lielab::domain::se & b, const size_t order){return Lielab::functions::dexpinv_numerical(a, b, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::so & a, const size_t order){return Lielab::functions::dexpinv_numerical(a, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::so & a, const Lielab::domain::so & b, const size_t order){return Lielab::functions::dexpinv_numerical(a, b, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::sp & a, const size_t order){return Lielab::functions::dexpinv_numerical(a, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::sp & a, const Lielab::domain::sp & b, const size_t order){return Lielab::functions::dexpinv_numerical(a, b, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::su & a, const size_t order){return Lielab::functions::dexpinv_numerical(a, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv_numerical", [](const Lielab::domain::su & a, const Lielab::domain::su & b, const size_t order){return Lielab::functions::dexpinv_numerical(a, b, order);}, "The numerical dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::cn & a, const size_t order){return Lielab::functions::dexpinv(a, order);}, "The dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::cn & a, const Lielab::domain::cn & b, const size_t order){return Lielab::functions::dexpinv(a, b, order);}, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::gl & a, const size_t order){return Lielab::functions::dexpinv(a, order);}, "The dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::gl & a, const Lielab::domain::gl & b, const size_t order){return Lielab::functions::dexpinv(a, b, order);}, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::glc & a, const size_t order){return Lielab::functions::dexpinv(a, order);}, "The dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::glc & a, const Lielab::domain::glc & b, const size_t order){return Lielab::functions::dexpinv(a, b, order);}, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::rn & a, const size_t order){return Lielab::functions::dexpinv(a, order);}, "The dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::rn & a, const Lielab::domain::rn & b, const size_t order){return Lielab::functions::dexpinv(a, b, order);}, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::se & a, const size_t order){return Lielab::functions::dexpinv(a, order);}, "The dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::se & a, const Lielab::domain::se & b, const size_t order){return Lielab::functions::dexpinv(a, b, order);}, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::so & a, const size_t order){return Lielab::functions::dexpinv(a, order);}, "The dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::so & a, const Lielab::domain::so & b, const size_t order){return Lielab::functions::dexpinv(a, b, order);}, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::sp & a, const size_t order){return Lielab::functions::dexpinv(a, order);}, "The dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::sp & a, const Lielab::domain::sp & b, const size_t order){return Lielab::functions::dexpinv(a, b, order);}, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::su & a, const size_t order){return Lielab::functions::dexpinv(a, order);}, "The dexpinv function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dexpinv", [](const Lielab::domain::su & a, const Lielab::domain::su & b, const size_t order){return Lielab::functions::dexpinv(a, b, order);}, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::cn & a, const size_t order){return Lielab::functions::dlog_numerical(a, order);}, "The numerical dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::cn & a, const Lielab::domain::cn & b, const size_t order){return Lielab::functions::dlog_numerical(a, b, order);}, "The numerical dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::gl & a, const size_t order){return Lielab::functions::dlog_numerical(a, order);}, "The numerical dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::gl & a, const Lielab::domain::gl & b, const size_t order){return Lielab::functions::dlog_numerical(a, b, order);}, "The numerical dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::glc & a, const size_t order){return Lielab::functions::dlog_numerical(a, order);}, "The numerical dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::glc & a, const Lielab::domain::glc & b, const size_t order){return Lielab::functions::dlog_numerical(a, b, order);}, "The numerical dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::rn & a, const size_t order){return Lielab::functions::dlog_numerical(a, order);}, "The numerical dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::rn & a, const Lielab::domain::rn & b, const size_t order){return Lielab::functions::dlog_numerical(a, b, order);}, "The numerical dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::se & a, const size_t order){return Lielab::functions::dlog_numerical(a, order);}, "The numerical dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::se & a, const Lielab::domain::se & b, const size_t order){return Lielab::functions::dlog_numerical(a, b, order);}, "The numerical dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::so & a, const size_t order){return Lielab::functions::dlog_numerical(a, order);}, "The numerical dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::so & a, const Lielab::domain::so & b, const size_t order){return Lielab::functions::dlog_numerical(a, b, order);}, "The numerical dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::sp & a, const size_t order){return Lielab::functions::dlog_numerical(a, order);}, "The numerical dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::sp & a, const Lielab::domain::sp & b, const size_t order){return Lielab::functions::dlog_numerical(a, b, order);}, "The numerical dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::su & a, const size_t order){return Lielab::functions::dlog_numerical(a, order);}, "The numerical dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog_numerical", [](const Lielab::domain::su & a, const Lielab::domain::su & b, const size_t order){return Lielab::functions::dlog_numerical(a, b, order);}, "The numerical dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::cn & a, const size_t order){return Lielab::functions::dlog(a, order);}, "The dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::cn & a, const Lielab::domain::cn & b, const size_t order){return Lielab::functions::dlog(a, b, order);}, "The dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::gl & a, const size_t order){return Lielab::functions::dlog(a, order);}, "The dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::gl & a, const Lielab::domain::gl & b, const size_t order){return Lielab::functions::dlog(a, b, order);}, "The dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::glc & a, const size_t order){return Lielab::functions::dlog(a, order);}, "The dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::glc & a, const Lielab::domain::glc & b, const size_t order){return Lielab::functions::dlog(a, b, order);}, "The dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::rn & a, const size_t order){return Lielab::functions::dlog(a, order);}, "The dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::rn & a, const Lielab::domain::rn & b, const size_t order){return Lielab::functions::dlog(a, b, order);}, "The dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::se & a, const size_t order){return Lielab::functions::dlog(a, order);}, "The dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::se & a, const Lielab::domain::se & b, const size_t order){return Lielab::functions::dlog(a, b, order);}, "The dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::so & a, const size_t order){return Lielab::functions::dlog(a, order);}, "The dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::so & a, const Lielab::domain::so & b, const size_t order){return Lielab::functions::dlog(a, b, order);}, "The dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::sp & a, const size_t order){return Lielab::functions::dlog(a, order);}, "The dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::sp & a, const Lielab::domain::sp & b, const size_t order){return Lielab::functions::dlog(a, b, order);}, "The dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::su & a, const size_t order){return Lielab::functions::dlog(a, order);}, "The dlog function.", py::arg("a"), py::arg("order") = 5);
    m_functions.def("dlog", [](const Lielab::domain::su & a, const Lielab::domain::su & b, const size_t order){return Lielab::functions::dlog(a, b, order);}, "The dlog function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("left_product", &Lielab::functions::left_product, "Default action by left product.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::cn>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::gl>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::glc>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::rn>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::se>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::so>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::sp>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::su>, "The pair function.");

    /*!
    * Begin content for the "kinematics" submodule.
    */
    py::module m_kinematics = m.def_submodule("kinematics", "The kinematics submodule.");
    m_kinematics.def("quaternions", &Lielab::kinematics::quaternions, "The quaternion kinematic equation of motion."); // TODO: Change this to use the actual quaternion from domain
    m_kinematics.def("eanglespace123", &Lielab::kinematics::eanglespace123, "The eanglespace123 equation of motion.");
    m_kinematics.def("eanglespace231", &Lielab::kinematics::eanglespace231, "The eanglespace231 equation of motion.");
    m_kinematics.def("eanglespace312", &Lielab::kinematics::eanglespace312, "The eanglespace312 equation of motion.");
    m_kinematics.def("eanglespace132", &Lielab::kinematics::eanglespace132, "The eanglespace132 equation of motion.");
    m_kinematics.def("eanglespace213", &Lielab::kinematics::eanglespace213, "The eanglespace213 equation of motion.");
    m_kinematics.def("eanglespace321", &Lielab::kinematics::eanglespace321, "The eanglespace321 equation of motion.");
    m_kinematics.def("eanglespace121", &Lielab::kinematics::eanglespace121, "The eanglespace121 equation of motion.");
    m_kinematics.def("eanglespace131", &Lielab::kinematics::eanglespace131, "The eanglespace131 equation of motion.");
    m_kinematics.def("eanglespace212", &Lielab::kinematics::eanglespace212, "The eanglespace212 equation of motion.");
    m_kinematics.def("eanglespace232", &Lielab::kinematics::eanglespace232, "The eanglespace232 equation of motion.");
    m_kinematics.def("eanglespace313", &Lielab::kinematics::eanglespace313, "The eanglespace313 equation of motion.");
    m_kinematics.def("eanglespace323", &Lielab::kinematics::eanglespace323, "The eanglespace323 equation of motion.");
    m_kinematics.def("eanglebody123", &Lielab::kinematics::eanglebody123, "The eanglebody123 equation of motion.");
    m_kinematics.def("eanglebody231", &Lielab::kinematics::eanglebody231, "The eanglebody231 equation of motion.");
    m_kinematics.def("eanglebody312", &Lielab::kinematics::eanglebody312, "The eanglebody312 equation of motion.");
    m_kinematics.def("eanglebody132", &Lielab::kinematics::eanglebody132, "The eanglebody132 equation of motion.");
    m_kinematics.def("eanglebody213", &Lielab::kinematics::eanglebody213, "The eanglebody213 equation of motion.");
    m_kinematics.def("eanglebody321", &Lielab::kinematics::eanglebody321, "The eanglebody321 equation of motion.");
    m_kinematics.def("eanglebody121", &Lielab::kinematics::eanglebody121, "The eanglebody121 equation of motion.");
    m_kinematics.def("eanglebody131", &Lielab::kinematics::eanglebody131, "The eanglebody131 equation of motion.");
    m_kinematics.def("eanglebody212", &Lielab::kinematics::eanglebody212, "The eanglebody212 equation of motion.");
    m_kinematics.def("eanglebody232", &Lielab::kinematics::eanglebody232, "The eanglebody232 equation of motion.");
    m_kinematics.def("eanglebody313", &Lielab::kinematics::eanglebody313, "The eanglebody313 equation of motion.");
    m_kinematics.def("eanglebody323", &Lielab::kinematics::eanglebody323, "The eanglebody323 equation of motion.");

    /*!
    * Begin content for the "optim" submodule.
    */
    py::module m_optim = m.def_submodule("optim", "The optim submodule.");

    py::class_<Lielab::optim::opt_golden>(m_optim, "opt_golden")
        .def(py::init())
        .def("init", &Lielab::optim::opt_golden::init)
        .def("step", &Lielab::optim::opt_golden::step)
        .def_readwrite("iterations", &Lielab::optim::opt_golden::iterations)
        .def_readwrite("max_iterations", &Lielab::optim::opt_golden::max_iterations)
        .def_readwrite("num_objective_evals", &Lielab::optim::opt_golden::num_objective_evals)
        .def_readwrite("num_jacobian_evals", &Lielab::optim::opt_golden::num_jacobian_evals)
        .def_readwrite("num_hessian_evals", &Lielab::optim::opt_golden::num_hessian_evals)
        .def_readwrite("algo_status", &Lielab::optim::opt_golden::algo_status)
        .def_readwrite("success", &Lielab::optim::opt_golden::success)
        .def_readwrite("tolerance", &Lielab::optim::opt_golden::tolerance)
        .def_readwrite("val_objective", &Lielab::optim::opt_golden::val_objective)
        .def_readwrite("val_jacobian", &Lielab::optim::opt_golden::val_jacobian)
        .def_readwrite("val_hessian", &Lielab::optim::opt_golden::val_hessian)
        .def_readwrite("x0", &Lielab::optim::opt_golden::x0)
        .def_readwrite("lower", &Lielab::optim::opt_golden::lower)
        .def_readwrite("upper", &Lielab::optim::opt_golden::upper)
        .def_readwrite("tau", &Lielab::optim::opt_golden::tau)
        .def_readwrite("_f1", &Lielab::optim::opt_golden::_f1)
        .def_readwrite("_f2", &Lielab::optim::opt_golden::_f2)
        .def_readwrite("_X", &Lielab::optim::opt_golden::_X)
        .def_readwrite("_X1", &Lielab::optim::opt_golden::_X1)
        .def_readwrite("_X2", &Lielab::optim::opt_golden::_X2)
        .def_readwrite("_A", &Lielab::optim::opt_golden::_A)
        .def_readwrite("_B", &Lielab::optim::opt_golden::_B);
    
    py::class_<Lielab::optim::hnewton>(m_optim, "hnewton")
        .def(py::init())
        .def("init", &Lielab::optim::hnewton::init)
        .def("step0", &Lielab::optim::hnewton::step0)
        .def("step1", &Lielab::optim::hnewton::step1)
        .def_readwrite("lower", &Lielab::optim::hnewton::lower)
        .def_readwrite("upper", &Lielab::optim::hnewton::upper)
        .def_readwrite("algo_status", &Lielab::optim::hnewton::algo_status)
        .def_readwrite("f", &Lielab::optim::hnewton::f)
        .def_readwrite("fnext", &Lielab::optim::hnewton::fnext)
        .def_readwrite("m", &Lielab::optim::hnewton::m)
        .def_readwrite("mnext", &Lielab::optim::hnewton::mnext);
    
    // py::class_<Lielab::optim::search_linearx>(m_optim, "search_linearx")
    //     .def(py::init())
    //     .def("init", &Lielab::optim::search_linearx::init)
    //     .def("step", &Lielab::optim::search_linearx::step)
    //     .def_readwrite("iterations", &Lielab::optim::search_linearx::iterations)
    //     .def_readwrite("max_iterations", &Lielab::optim::search_linearx::max_iterations)
    //     .def_readwrite("num_objective_evals", &Lielab::optim::search_linearx::num_objective_evals)
    //     .def_readwrite("num_jacobian_evals", &Lielab::optim::search_linearx::num_jacobian_evals)
    //     .def_readwrite("num_hessian_evals", &Lielab::optim::search_linearx::num_hessian_evals)
    //     .def_readwrite("algo_status", &Lielab::optim::search_linearx::algo_status)
    //     .def_readwrite("success", &Lielab::optim::search_linearx::success)
    //     .def_readwrite("tolerance", &Lielab::optim::search_linearx::tolerance)
    //     .def_readwrite("val_objective", &Lielab::optim::search_linearx::val_objective)
    //     .def_readwrite("val_jacobian", &Lielab::optim::search_linearx::val_jacobian)
    //     .def_readwrite("val_hessian", &Lielab::optim::search_linearx::val_hessian)
    //     .def_readwrite("_dx", &Lielab::optim::search_linearx::_dx)
    //     .def_readwrite("_x", &Lielab::optim::search_linearx::_x)
    //     .def_readwrite("_x1", &Lielab::optim::search_linearx::_x1)
    //     .def_readwrite("_x2", &Lielab::optim::search_linearx::_x2)
    //     .def_readwrite("_y1", &Lielab::optim::search_linearx::_y1)
    //     .def_readwrite("_y2", &Lielab::optim::search_linearx::_y2)
    //     .def_readwrite("k", &Lielab::optim::search_linearx::k)
    //     .def_readwrite("_lo", &Lielab::optim::search_linearx::_lo)
    //     .def_readwrite("_hi", &Lielab::optim::search_linearx::_hi)
    //     .def_readwrite("fdx", &Lielab::optim::search_linearx::fdx)
    //     .def_readwrite("lower", &Lielab::optim::search_linearx::lower)
    //     .def_readwrite("upper", &Lielab::optim::search_linearx::upper);
    

    /*!
    * Begin content for the "transform" submodule.
    */

    py::module m_transform = m.def_submodule("transform", "The transform submodule.");
    m_transform.def("dcm_to_quaternion", &Lielab::transform::dcm_to_quaternion, "The dcm_to_quaternion function.");
    m_transform.def("dcm_to_eanglebody123", &Lielab::transform::dcm_to_eanglebody123<double>, "The dcm_to_eanglebody123 function.");
    m_transform.def("dcm_to_eanglebody231", &Lielab::transform::dcm_to_eanglebody231<double>, "The dcm_to_eanglebody231 function.");
    m_transform.def("dcm_to_eanglebody312", &Lielab::transform::dcm_to_eanglebody312<double>, "The dcm_to_eanglebody312 function.");
    m_transform.def("dcm_to_eanglebody132", &Lielab::transform::dcm_to_eanglebody132<double>, "The dcm_to_eanglebody132 function.");
    m_transform.def("dcm_to_eanglebody213", &Lielab::transform::dcm_to_eanglebody213<double>, "The dcm_to_eanglebody213 function.");
    m_transform.def("dcm_to_eanglebody321", &Lielab::transform::dcm_to_eanglebody321<double>, "The dcm_to_eanglebody321 function.");
    m_transform.def("dcm_to_eanglebody121", &Lielab::transform::dcm_to_eanglebody121<double>, "The dcm_to_eanglebody121 function.");
    m_transform.def("dcm_to_eanglebody131", &Lielab::transform::dcm_to_eanglebody131<double>, "The dcm_to_eanglebody131 function.");
    m_transform.def("dcm_to_eanglebody212", &Lielab::transform::dcm_to_eanglebody212<double>, "The dcm_to_eanglebody212 function.");
    m_transform.def("dcm_to_eanglebody232", &Lielab::transform::dcm_to_eanglebody232<double>, "The dcm_to_eanglebody232 function.");
    m_transform.def("dcm_to_eanglebody313", &Lielab::transform::dcm_to_eanglebody313<double>, "The dcm_to_eanglebody313 function.");
    m_transform.def("dcm_to_eanglebody323", &Lielab::transform::dcm_to_eanglebody323<double>, "The dcm_to_eanglebody323 function.");
    m_transform.def("dcm_to_eanglespace123", &Lielab::transform::dcm_to_eanglespace123<double>, "The dcm_to_eanglespace123 function.");
    m_transform.def("dcm_to_eanglespace231", &Lielab::transform::dcm_to_eanglespace231<double>, "The dcm_to_eanglespace231 function.");
    m_transform.def("dcm_to_eanglespace312", &Lielab::transform::dcm_to_eanglespace312<double>, "The dcm_to_eanglespace312 function.");
    m_transform.def("dcm_to_eanglespace132", &Lielab::transform::dcm_to_eanglespace132<double>, "The dcm_to_eanglespace132 function.");
    m_transform.def("dcm_to_eanglespace213", &Lielab::transform::dcm_to_eanglespace213<double>, "The dcm_to_eanglespace213 function.");
    m_transform.def("dcm_to_eanglespace321", &Lielab::transform::dcm_to_eanglespace321<double>, "The dcm_to_eanglespace321 function.");
    m_transform.def("dcm_to_eanglespace121", &Lielab::transform::dcm_to_eanglespace121<double>, "The dcm_to_eanglespace121 function.");
    m_transform.def("dcm_to_eanglespace131", &Lielab::transform::dcm_to_eanglespace131<double>, "The dcm_to_eanglespace131 function.");
    m_transform.def("dcm_to_eanglespace212", &Lielab::transform::dcm_to_eanglespace212<double>, "The dcm_to_eanglespace212 function.");
    m_transform.def("dcm_to_eanglespace232", &Lielab::transform::dcm_to_eanglespace232<double>, "The dcm_to_eanglespace232 function.");
    m_transform.def("dcm_to_eanglespace313", &Lielab::transform::dcm_to_eanglespace313<double>, "The dcm_to_eanglespace313 function.");
    m_transform.def("dcm_to_eanglespace323", &Lielab::transform::dcm_to_eanglespace323<double>, "The dcm_to_eanglespace323 function.");
    m_transform.def("dcm_to_gibbs", &Lielab::transform::dcm_to_gibbs, "DCM to Gibbs.");
    m_transform.def("quaternion_to_dcm", &Lielab::transform::quaternion_to_dcm, "The quaternion_to_dcm function.");

    py::module m_topos = m.def_submodule("topos", "The topos submodule.");

    m_topos.def("Ad", &Lielab::topos::Ad);
    m_topos.def("cayley1", &Lielab::topos::cayley1);
    m_topos.def("dcayley1inv", &Lielab::topos::dcayley1inv);
    m_topos.def("dexp", &Lielab::topos::dexp, "The dexp function.");
    m_topos.def("dexpinv", &Lielab::topos::dexpinv);
    m_topos.def("exp", &Lielab::topos::exp);
    m_topos.def("log", &Lielab::topos::log);
    
    // py::enum_<Lielab::topos::RKTYPE>(m_topos, "RKTYPE")
    //     .value("RKTYPE_NONE", Lielab::topos::RKTYPE::RKTYPE_NONE)
    //     .value("RKTYPE_EXPLICIT", Lielab::topos::RKTYPE::RKTYPE_EXPLICIT)
    //     .value("RKTYPE_IMPLICIT", Lielab::topos::RKTYPE::RKTYPE_IMPLICIT);
    
    py::class_<Lielab::topos::IntegralCurve>(m_topos, "IntegralCurve")
        .def(py::init<>())
        .def(py::init<Lielab::topos::IntegralCurve>())
        .def(py::init<size_t>())
        .def_readwrite("chunk", &Lielab::topos::IntegralCurve::chunk)
        .def_readwrite("length", &Lielab::topos::IntegralCurve::length)
        .def_readwrite("num_eoms", &Lielab::topos::IntegralCurve::num_eoms)
        .def_readwrite("chs", &Lielab::topos::IntegralCurve::chs)
        .def_readwrite("t", &Lielab::topos::IntegralCurve::t)
        .def_readwrite("y", &Lielab::topos::IntegralCurve::y);
    
    py::class_<Lielab::topos::TSOutput>(m_topos, "TSOutput")
        .def(py::init())
        .def(py::init<Lielab::domain::CompositeManifold>())
        .def(py::init<Lielab::domain::CompositeManifold, Lielab::domain::CompositeManifold, double>())
        .def_readonly("low", &Lielab::topos::TSOutput::low)
        .def_readonly("high", &Lielab::topos::TSOutput::high)
        .def_readonly("error", &Lielab::topos::TSOutput::error);

    py::enum_<Lielab::topos::RKMETHOD>(m_topos, "RKMETHOD")
        .value("E1", Lielab::topos::RKMETHOD::E1)
        .value("RK45", Lielab::topos::RKMETHOD::RK45);
    
    py::class_<Lielab::topos::TimeStepper>(m_topos, "TimeStepper")
        .def_readwrite("algo_status", &Lielab::topos::TimeStepper::algo_status)
        .def_readwrite("iterations", &Lielab::topos::TimeStepper::iterations)
        .def_readwrite("max_iterations", &Lielab::topos::TimeStepper::max_iterations)
        .def_readwrite("success", &Lielab::topos::TimeStepper::success)
        .def_readonly("A", &Lielab::topos::TimeStepper::A)
        .def_readonly("B", &Lielab::topos::TimeStepper::B)
        .def_readonly("Bhat", &Lielab::topos::TimeStepper::Bhat)
        .def_readonly("C", &Lielab::topos::TimeStepper::C)
        .def_readonly("n", &Lielab::topos::TimeStepper::n)
        .def_readonly("order", &Lielab::topos::TimeStepper::order)
        .def_readonly("variable_step", &Lielab::topos::TimeStepper::variable_step)
        .def(py::init())
        .def(py::init<Lielab::topos::RKMETHOD>())
        .def("init", &Lielab::topos::TimeStepper::init)
        .def("step", &Lielab::topos::TimeStepper::step);

    py::class_<Lielab::topos::MuntheKaas>(m_topos, "MuntheKaas")
        .def_readwrite("algo_status", &Lielab::topos::MuntheKaas::algo_status)
        .def_readwrite("iterations", &Lielab::topos::TimeStepper::iterations)
        .def_readwrite("max_iterations", &Lielab::topos::TimeStepper::max_iterations)
        .def_readwrite("success", &Lielab::topos::TimeStepper::success)
        .def_readonly("_KK", &Lielab::topos::MuntheKaas::_KK)
        .def_readonly("_U", &Lielab::topos::MuntheKaas::_U)
        .def_readonly("A", &Lielab::topos::MuntheKaas::A)
        .def_readonly("B", &Lielab::topos::MuntheKaas::B)
        .def_readonly("Bhat", &Lielab::topos::MuntheKaas::Bhat)
        .def_readonly("C", &Lielab::topos::MuntheKaas::C)
        .def_readonly("n", &Lielab::topos::MuntheKaas::n)
        .def_readonly("order", &Lielab::topos::MuntheKaas::order)
        .def_readonly("variable_step", &Lielab::topos::MuntheKaas::variable_step)
        .def(py::init())
        .def(py::init<Lielab::topos::RKMETHOD>())
        .def("init", &Lielab::topos::MuntheKaas::init)
        .def("step_0", &Lielab::topos::MuntheKaas::step_0)
        .def("step_1", &Lielab::topos::MuntheKaas::step_1)
        .def("set_dy", &Lielab::topos::MuntheKaas::set_dy)
        .def("postprocess", &Lielab::topos::MuntheKaas::postprocess)
        .def_readwrite("next_t", &Lielab::topos::MuntheKaas::next_t)
        .def_readwrite("next_y", &Lielab::topos::MuntheKaas::next_y)
        .def("__repr__",
        [](const Lielab::topos::MuntheKaas & self)
        {
            return "<lielab.topos.MuntheKaas>";
        })
        .def("__str__",
        [](const Lielab::topos::MuntheKaas & self)
        {
            return "<lielab.topos.MuntheKaas>";
        });
    
    py::class_<Lielab::topos::Flow>(m_topos, "Flow")
        .def(py::init())
        .def("init", &Lielab::topos::Flow::init)
        .def("step0", &Lielab::topos::Flow::step0)
        .def("step", &Lielab::topos::Flow::step)
        .def("stepE", &Lielab::topos::Flow::stepE)
        .def("postprocess", &Lielab::topos::Flow::postprocess)
        .def("new_step_size", &Lielab::topos::Flow::new_step_size)
        .def_readwrite("algo_status", &Lielab::topos::Flow::algo_status)
        .def_readwrite("_ynext", &Lielab::topos::Flow::_ynext)
        .def_readwrite("iterations", &Lielab::topos::Flow::iterations)
        .def_readwrite("small", &Lielab::topos::Flow::small)
        .def_readwrite("large", &Lielab::topos::Flow::large)
        .def_readwrite("pessimist", &Lielab::topos::Flow::pessimist)
        .def_readwrite("accept", &Lielab::topos::Flow::accept)
        .def_readwrite("tol", &Lielab::topos::Flow::tol)
        .def_readwrite("dt_min", &Lielab::topos::Flow::dt_min)
        .def_readwrite("dt_max", &Lielab::topos::Flow::dt_max)
        .def_readwrite("dt", &Lielab::topos::Flow::dt)
        .def_readwrite("default_local", &Lielab::topos::Flow::default_local)
        .def_readwrite("default_global", &Lielab::topos::Flow::default_global)
        .def_readwrite("variable_time_step", &Lielab::topos::Flow::variable_time_step)
    //     .def_readwrite("new_exact", &Lielab::topos::Flow::new_exact)
    //     .def_readwrite("num_step", &Lielab::topos::Flow::num_step)
        .def_readwrite("_out", &Lielab::topos::Flow::_out)
        .def_readwrite("_dt", &Lielab::topos::Flow::_dt)
        .def_readwrite("stepper", &Lielab::topos::Flow::stepper)
        .def("__repr__",
        [](const Lielab::topos::Flow & self)
        {
            return "<lielab.topos.Flow>";
        })
        .def("__str__",
        [](const Lielab::topos::Flow & self)
        {
            return "<lielab.topos.Flow>";
        });
}
