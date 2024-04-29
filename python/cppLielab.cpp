#include <lielab>
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
    m.attr("__author__") = lielab::AUTHOR;
    m.attr("__contact__") = lielab::CONTACT;
    m.attr("__location__") = lielab::LOCATION;
    m.attr("__version__") = lielab::VERSION;

    /*!
    * Begin content for the "lielab" module.
    */

    py::enum_<lielab::ALGO_STATUS>(m, "ALGO_STATUS")
        .value("OK", lielab::ALGO_STATUS::OK)
        .value("MAXITER", lielab::ALGO_STATUS::MAXITER)
        .value("FINISHED", lielab::ALGO_STATUS::FINISHED);

    /*!
    * Begin content for the "domain" submodule.
    */
    py::module m_domain = m.def_submodule("domain", "The domain submodule.");

    /*!
    * Bindings for lielab::domain::gl
    */

    py::class_<lielab::domain::gl>(m_domain, "gl")
        .def_readwrite("_data", &lielab::domain::gl::_data)
        .def_readonly_static("abelian", &lielab::domain::gl::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXd &>())
        .def("basis", &lielab::domain::gl::basis)
        .def("project", &lielab::domain::gl::project)
        .def("get_dimension", &lielab::domain::gl::get_dimension)
        .def("get_vector", &lielab::domain::gl::get_vector)
        .def("get_ados_representation", &lielab::domain::gl::get_ados_representation)
        .def("set_vector", &lielab::domain::gl::set_vector)
        .def("__call__", [](lielab::domain::gl & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const lielab::domain::gl & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(-py::self)
        .def(py::self * int())
        .def(py::self * double())
        .def(py::self * py::self)
        .def(py::self *= int())
        .def(py::self *= double())
        .def(py::self *= py::self)
        .def(py::self / int())
        .def(py::self / double())
        .def(py::self /= int())
        .def(py::self /= double())
        .def(int() * py::self)
        .def(double() * py::self)
        .def("__repr__", [](const lielab::domain::gl & self)
        {
            return "<lielab.domain.gl>";
        })
        .def("__str__", [](const lielab::domain::gl & self)
        {
            return matstr(self._data);
        });

    /*!
    * Bindings for lielab::domain::rn
    */

    py::class_<lielab::domain::rn>(m_domain, "rn")
        .def_readwrite("_data", &lielab::domain::rn::_data)
        .def_readwrite("shape", &lielab::domain::rn::shape)
        .def_readonly_static("abelian", &lielab::domain::rn::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXd &>())
        .def("basis", &lielab::domain::rn::basis)
        .def("project", &lielab::domain::rn::project)
        .def("get_dimension", &lielab::domain::rn::get_dimension)
        .def("get_vector", &lielab::domain::rn::get_vector)
        .def("get_ados_representation", &lielab::domain::rn::get_ados_representation)
        .def("set_vector", &lielab::domain::rn::set_vector)
        .def("__call__", [](const lielab::domain::rn & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const lielab::domain::rn & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(-py::self)
        .def(py::self * int())
        .def(py::self * double())
        .def(py::self * py::self)
        .def(py::self *= int())
        .def(py::self *= double())
        .def(py::self *= py::self)
        .def(py::self / int())
        .def(py::self / double())
        .def(py::self /= int())
        .def(py::self /= double())
        .def(int() * py::self)
        .def(double() * py::self)
        .def("__repr__", [](const lielab::domain::rn & self)
        {
            return "<lielab.domain.rn>";
        })
        .def("__str__", [](const lielab::domain::rn & self)
        {
            return matstr(self._data);
        });


    /*!
    * Bindings for lielab::domain::so
    */

    py::class_<lielab::domain::so>(m_domain, "so")
        .def_readwrite("_data", &lielab::domain::so::_data)
        .def_readwrite("shape", &lielab::domain::so::shape)
        .def_readonly_static("abelian", &lielab::domain::so::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXd &>())
        .def("basis", &lielab::domain::so::basis)
        .def("project", &lielab::domain::so::project)
        .def("get_dimension", &lielab::domain::so::get_dimension)
        .def("get_vector", &lielab::domain::so::get_vector)
        .def("get_ados_representation", &lielab::domain::so::get_ados_representation)
        .def("set_vector", &lielab::domain::so::set_vector)
        .def("__call__", [](const lielab::domain::so & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const lielab::domain::so & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(-py::self)
        .def(py::self * int())
        .def(py::self * double())
        .def(py::self * py::self)
        .def(py::self *= int())
        .def(py::self *= double())
        .def(py::self *= py::self)
        .def(py::self / int())
        .def(py::self / double())
        .def(py::self /= int())
        .def(py::self /= double())
        .def(int() * py::self)
        .def(double() * py::self)
        .def("__repr__", [](const lielab::domain::so & self)
        {
            return "<lielab.domain.so>";
        })
        .def("__str__", [](const lielab::domain::so & self)
        {
            return matstr(self._data);
        });

    
    /*!
    * Bindings for lielab::domain::sp
    */
    
    py::class_<lielab::domain::sp>(m_domain, "sp")
        .def_readwrite("_data", &lielab::domain::sp::_data)
        .def_readwrite("shape", &lielab::domain::sp::shape)
        .def_readonly_static("abelian", &lielab::domain::sp::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXd &>())
        .def("basis", &lielab::domain::sp::basis)
        .def("project", &lielab::domain::sp::project)
        .def("get_dimension", &lielab::domain::sp::get_dimension)
        .def("get_vector", &lielab::domain::sp::get_vector)
        .def("get_ados_representation", &lielab::domain::sp::get_ados_representation)
        .def("set_vector", &lielab::domain::sp::set_vector)
        .def("__call__", [](const lielab::domain::sp & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const lielab::domain::sp & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(-py::self)
        .def(py::self * int())
        .def(py::self * double())
        .def(py::self * py::self)
        .def(py::self *= int())
        .def(py::self *= double())
        .def(py::self *= py::self)
        .def(py::self / int())
        .def(py::self / double())
        .def(py::self /= int())
        .def(py::self /= double())
        .def(int() * py::self)
        .def(double() * py::self)
        .def("__repr__", [](const lielab::domain::sp & self)
        {
            return "<lielab.domain.sp>";
        })
        .def("__str__", [](const lielab::domain::sp & self)
        {
            return matstr(self._data);
        });
    
    
    /*!
    * Bindings for lielab::domain::su
    */
    
    py::class_<lielab::domain::su>(m_domain, "su")
        .def_readwrite("_data", &lielab::domain::su::_data)
        .def_readwrite("shape", &lielab::domain::su::shape)
        .def_readonly_static("abelian", &lielab::domain::su::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXcd &>())
        .def("basis", &lielab::domain::su::basis)
        // .def("project", &lielab::domain::su::project) // TODO: Implement project
        .def("get_dimension", &lielab::domain::su::get_dimension)
        .def("get_vector", &lielab::domain::su::get_vector)
        .def("get_ados_representation", &lielab::domain::su::get_ados_representation)
        .def("set_vector", &lielab::domain::su::set_vector)
        .def("__call__", [](const lielab::domain::su & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const lielab::domain::su & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(-py::self)
        .def(py::self * int())
        .def(py::self * double())
        // .def(py::self * std::complex<int>()) // TODO: Complex integers seem bugged in Eigen right now
        .def(py::self * std::complex<double>())
        .def(py::self * py::self)
        .def(py::self *= int())
        .def(py::self *= double())
        // .def(py::self *= std::complex<int>()) // TODO: Complex integers seem bugged in Eigen right now
        .def(py::self *= std::complex<double>())
        .def(py::self *= py::self)
        .def(py::self / int())
        .def(py::self / double())
        // .def(py::self / std::complex<int>())
        .def(py::self / std::complex<double>())
        .def(py::self /= int())
        .def(py::self /= double())
        // .def(py::self /= std::complex<int>())
        .def(py::self /= std::complex<double>())
        .def(int() * py::self)
        .def(double() * py::self)
        // .def(std::complex<int>() * py::self)
        .def(std::complex<double>() * py::self)
        .def("__repr__", [](const lielab::domain::su & self)
        {
            return "<lielab.domain.su>";
        })
        .def("__str__", [](const lielab::domain::su & self)
        {
            return matstr(self._data);
        });
    
    /*!
     * Lie Groups
     */
    
    py::class_<lielab::domain::GL>(m_domain, "GL")
        .def_readwrite("_data", &lielab::domain::GL::_data)
        .def_readonly_static("abelian", &lielab::domain::GL::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXd>())
        .def("project", &lielab::domain::GL::project)
        .def("get_dimension", &lielab::domain::GL::get_dimension)
        .def("get_ados_representation", &lielab::domain::GL::get_ados_representation)
        .def("inverse", &lielab::domain::GL::inverse)
        .def("serialize", &lielab::domain::GL::serialize)
        .def("unserialize", &lielab::domain::GL::unserialize)
        .def("__call__", [](const lielab::domain::GL & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const lielab::domain::GL & self)
        {
            return "<lielab.domain.GL>";
        })
        .def("__str__", [](const lielab::domain::GL & self)
        {
            return matstr(self._data);
        });

    py::class_<lielab::domain::RN>(m_domain, "RN")
        .def_readwrite("_data", &lielab::domain::RN::_data)
        .def_readwrite("shape", &lielab::domain::RN::shape)
        .def_readonly_static("abelian", &lielab::domain::RN::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXd>())
        .def("project", &lielab::domain::RN::project)
        .def("get_dimension", &lielab::domain::RN::get_dimension)
        .def("get_ados_representation", &lielab::domain::RN::get_ados_representation)
        .def("inverse", &lielab::domain::RN::inverse)
        .def("serialize", &lielab::domain::RN::serialize)
        .def("unserialize", &lielab::domain::RN::unserialize)
        .def("__call__", [](const lielab::domain::RN & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const lielab::domain::RN & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const lielab::domain::RN & self)
        {
            return "<lielab.domain.RN>";
        })
        .def("__str__", [](const lielab::domain::RN & self)
        {
            return matstr(self._data);
        });


    py::class_<lielab::domain::SO>(m_domain, "SO")
        .def_readwrite("_data", &lielab::domain::SO::_data)
        .def_readwrite("shape", &lielab::domain::SO::shape)
        .def_readonly_static("abelian", &lielab::domain::SO::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXd>())
        .def("project", &lielab::domain::SO::project)
        .def("get_dimension", &lielab::domain::SO::get_dimension)
        .def("get_ados_representation", &lielab::domain::SO::get_ados_representation)
        .def("inverse", &lielab::domain::SO::inverse)
        .def("serialize", &lielab::domain::SO::serialize)
        .def("unserialize", &lielab::domain::SO::unserialize)
        .def("__call__", [](const lielab::domain::SO & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const lielab::domain::SO & self)
        {
            return "<lielab.domain.SO>";
        })
        .def("__str__", [](const lielab::domain::SO & self)
        {
            return matstr(self._data);
        });
    

    py::class_<lielab::domain::SP>(m_domain, "SP")
        .def_readwrite("_data", &lielab::domain::SP::_data)
        .def_readwrite("shape", &lielab::domain::SP::shape)
        .def_readonly_static("abelian", &lielab::domain::SP::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXd>())
        // .def("project", &lielab::domain::SP::project) // TODO:
        .def("get_dimension", &lielab::domain::SP::get_dimension)
        .def("get_ados_representation", &lielab::domain::SP::get_ados_representation)
        .def("inverse", &lielab::domain::SP::inverse)
        .def("serialize", &lielab::domain::SP::serialize)
        // .def("unserialize", &lielab::domain::SP::unserialize) // TODO:
        .def("__call__", [](const lielab::domain::SP & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const lielab::domain::SP & self)
        {
            return "<lielab.domain.SP>";
        })
        .def("__str__", [](const lielab::domain::SP & self)
        {
            return matstr(self._data);
        });

    py::class_<lielab::domain::SU>(m_domain, "SU")
        .def_readwrite("_data", &lielab::domain::SU::_data)
        .def_readwrite("shape", &lielab::domain::SU::shape)
        .def_readonly_static("abelian", &lielab::domain::SU::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXcd>())
        .def_static("Quaternion", py::overload_cast<>(&lielab::domain::SU::Quaternion))
        .def_static("Quaternion", py::overload_cast<const double, const double, const double, const double>(&lielab::domain::SU::Quaternion))
        // .def("project", &lielab::domain::SU::project) // TODO:
        // .def("get_dimension", &lielab::domain::SU::get_dimension) // TODO:
        .def("get_ados_representation", &lielab::domain::SU::get_ados_representation)
        .def("inverse", &lielab::domain::SU::inverse)
        .def("serialize", &lielab::domain::SU::serialize)
        .def("unserialize", &lielab::domain::SU::unserialize)
        .def("__call__", [](const lielab::domain::SU & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const lielab::domain::SU & self)
        {
            return "<lielab.domain.SU>";
        })
        .def("__str__", [](const lielab::domain::SU & self)
        {
            return matstr(self._data);
        });

    py::class_<lielab::domain::halie>(m_domain, "halie")
        .def(py::init())
        .def(py::init<std::vector<lielab::domain::halie::TYPES>>())
        .def_readwrite("space", &lielab::domain::halie::space)
        .def("get_dimension", &lielab::domain::halie::get_dimension)
        .def("get_shape", &lielab::domain::halie::get_shape)
        .def("get_vector", &lielab::domain::halie::get_vector)
        .def("set_vector", &lielab::domain::halie::set_vector)
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(-py::self)
        .def(py::self * int())
        .def(py::self * double())
        .def(py::self * py::self)
        .def(py::self *= int())
        .def(py::self *= double())
        .def(py::self *= py::self)
        .def(py::self / int())
        .def(py::self / double())
        .def(py::self /= int())
        .def(py::self /= double())
        .def(int() * py::self)
        .def(double() * py::self)
        .def("__repr__",
        [](const lielab::domain::halie & self)
        {
            return "<lielab.domain.halie>";
        })
        .def("__str__",
        [](const lielab::domain::halie & self)
        {
            return self.to_string();
        });
    

    py::class_<lielab::domain::hmlie>(m_domain, "hmlie")
        .def(py::init())
        .def(py::init<std::vector<lielab::domain::hmlie::TYPES>>())
        .def_readwrite("space", &lielab::domain::hmlie::space)
        .def("get_shape", &lielab::domain::hmlie::get_shape)
        .def("serialize", &lielab::domain::hmlie::serialize)
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("inverse", &lielab::domain::hmlie::inverse)
        .def("__repr__",
        [](const lielab::domain::hmlie & self)
        {
            return "<lielab.domain.hmlie>";
        })
        .def("__str__",
        [](const lielab::domain::hmlie & self)
        {
            return self.to_string();
        });

    
    py::module m_dynamics = m.def_submodule("dynamics", "The dynamics submodule.");
    m_dynamics.def("vfex1", &lielab::dynamics::vfex1);
    m_dynamics.def("vfex2", &lielab::dynamics::vfex2);
    m_dynamics.def("newton_eom", &lielab::dynamics::newton_eom);
    m_dynamics.def("newton_action", &lielab::dynamics::newton_action);


    /*!
    * Begin content for the "functions" submodule.
    */

    py::module m_functions = m.def_submodule("functions", "The functions submodule.");
    // m_functions.def("pair", &lielab::functions::pair<lielab::domain::gl>, "The pair function.");
    m_functions.def("pair", &lielab::functions::pair<lielab::domain::rn>, "The pair function.");
    m_functions.def("pair", &lielab::functions::pair<lielab::domain::so>, "The pair function.");
    m_functions.def("pair", &lielab::functions::pair<lielab::domain::sp>, "The pair function.");
    m_functions.def("pair", &lielab::functions::pair<lielab::domain::su>, "The pair function.");
    // m_functions.def("copair", &lielab::functions::copair<lielab::domain::gl>, "The copair function.");
    // m_functions.def("copair", &lielab::functions::copair<lielab::domain::rn>, "The copair function.");
    // m_functions.def("copair", &lielab::functions::copair<lielab::domain::so>, "The copair function.");
    // m_functions.def("copair", &lielab::functions::copair<lielab::domain::sp>, "The copair function.");
    // m_functions.def("copair", &lielab::functions::copair<lielab::domain::su>, "The copair function.");
    m_functions.def("factorial", &lielab::functions::factorial, "The factorial function");
    m_functions.def("Ad", &lielab::functions::Ad<lielab::domain::gl>, "The Ad function.");
    m_functions.def("Ad", &lielab::functions::Ad<lielab::domain::rn>, "The Ad function.");
    m_functions.def("Ad", &lielab::functions::Ad<lielab::domain::so>, "The Ad function.");
    m_functions.def("Ad", &lielab::functions::Ad<lielab::domain::sp>, "The Ad function.");
    m_functions.def("Ad", &lielab::functions::Ad<lielab::domain::su>, "The Ad function.");
    m_functions.def("commutator", &lielab::functions::commutator<lielab::domain::gl>, "The commutator function.");
    m_functions.def("commutator", &lielab::functions::commutator<lielab::domain::rn>, "The commutator function.");
    m_functions.def("commutator", &lielab::functions::commutator<lielab::domain::so>, "The commutator function.");
    m_functions.def("commutator", &lielab::functions::commutator<lielab::domain::sp>, "The commutator function.");
    m_functions.def("commutator", &lielab::functions::commutator<lielab::domain::su>, "The commutator function.");
    // m_functions.def("cayley1", &lielab::functions::cayley1<lielab::domain::gl>, "The cayley1 function.");
    m_functions.def("cayley1", &lielab::functions::cayley1<lielab::domain::rn>, "The cayley1 function.");
    m_functions.def("cayley1", &lielab::functions::cayley1<lielab::domain::so>, "The cayley1 function.");
    m_functions.def("cayley1", &lielab::functions::cayley1<lielab::domain::sp>, "The cayley1 function.");
    m_functions.def("cayley1", &lielab::functions::cayley1<lielab::domain::su>, "The cayley1 function.");
    // m_functions.def("cayley1", &lielab::functions::cayley2<lielab::domain::gl>, "The cayley2 function.");
    m_functions.def("cayley2", &lielab::functions::cayley2<lielab::domain::rn>, "The cayley2 function.");
    m_functions.def("cayley2", &lielab::functions::cayley2<lielab::domain::so>, "The cayley2 function.");
    m_functions.def("cayley2", &lielab::functions::cayley2<lielab::domain::sp>, "The cayley2 function.");
    // m_functions.def("cayley2", &lielab::functions::cayley2<lielab::domain::su>, "The cayley2 function.");
    m_functions.def("Killing", &lielab::functions::Killing<lielab::domain::rn>, "The Killing function.");
    m_functions.def("Killing", &lielab::functions::Killing<lielab::domain::so>, "The Killing function.");
    m_functions.def("Killing", &lielab::functions::Killing<lielab::domain::sp>, "The Killing function.");
    // m_functions.def("Killing", &lielab::functions::Killing<lielab::domain::su>, "The Killing function.");
    m_functions.def("Killingform", &lielab::functions::Killingform<lielab::domain::gl>, "The Killingform function.");
    m_functions.def("Killingform", &lielab::functions::Killingform<lielab::domain::rn>, "The Killingform function.");
    m_functions.def("Killingform", &lielab::functions::Killingform<lielab::domain::so>, "The Killingform function.");
    m_functions.def("Killingform", &lielab::functions::Killingform<lielab::domain::sp>, "The Killingform function.");
    // m_functions.def("Killingform", &lielab::functions::Killingform<lielab::domain::su>, "The Killingform function.");
    // m_functions.def("ad", &lielab::functions::ad<lielab::domain::gl>, "The ad function.");
    m_functions.def("ad", &lielab::functions::ad<lielab::domain::rn>, "The ad function.");
    m_functions.def("ad", &lielab::functions::ad<lielab::domain::so>, "The ad function.");
    m_functions.def("ad", &lielab::functions::ad<lielab::domain::sp>, "The ad function.");
    m_functions.def("ad", &lielab::functions::ad<lielab::domain::su>, "The ad function.");
    // m_functions.def("coAd", &lielab::functions::coAd<lielab::domain::GL>, "The coAd function.");
    // m_functions.def("coAd", &lielab::functions::coAd<lielab::domain::RN>, "The coAd function.");
    // m_functions.def("coAd", &lielab::functions::coAd<lielab::domain::SO>, "The coAd function.");
    // m_functions.def("coAd", &lielab::functions::coAd<lielab::domain::SP>, "The coAd function.");
    // m_functions.def("coAd", &lielab::functions::coAd<lielab::domain::SU>, "The coAd function."); // SU needs get_dimension()
    // m_functions.def("coad", &lielab::functions::coad<lielab::domain::gl>, "The coad function.");
    // m_functions.def("coad", &lielab::functions::coad<lielab::domain::rn>, "The coad function.");
    // m_functions.def("coad", &lielab::functions::coad<lielab::domain::so>, "The coad function.");
    // m_functions.def("coad", &lielab::functions::coad<lielab::domain::sp>, "The coad function.");
    // m_functions.def("coad", &lielab::functions::coad<lielab::domain::su>, "The coad function."); // TODO: su needs basis()
    m_functions.def("exp", &lielab::functions::exp<lielab::domain::gl>, "The exponential function.");
    m_functions.def("exp", &lielab::functions::exp<lielab::domain::rn>, "The exponential function.");
    m_functions.def("exp", &lielab::functions::exp<lielab::domain::so>, "The exponential function.");
    m_functions.def("exp", &lielab::functions::exp<lielab::domain::sp>, "The exponential function.");
    m_functions.def("exp", &lielab::functions::exp<lielab::domain::su>, "The exponential function.");
    m_functions.def("log", &lielab::functions::log<lielab::domain::GL>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &lielab::functions::log<lielab::domain::RN>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &lielab::functions::log<lielab::domain::SO>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &lielab::functions::log<lielab::domain::SP>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &lielab::functions::log<lielab::domain::SU>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("bernoulli", &lielab::functions::bernoulli, "The bernoulli function.");
    // m_functions.def("log", &lielab::functions::log<lielab::domain::Quaternion>, "The logarithm function."); // TODO: Correctly implement this
    m_functions.def("dcayley1inv", &lielab::functions::dcayley1inv<lielab::domain::gl>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &lielab::functions::dcayley1inv<lielab::domain::rn>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &lielab::functions::dcayley1inv<lielab::domain::so>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &lielab::functions::dcayley1inv<lielab::domain::sp>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &lielab::functions::dcayley1inv<lielab::domain::su>, "The dcayley1inv function.");
    m_functions.def("dexp", &lielab::functions::dexp<lielab::domain::gl>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", &lielab::functions::dexp<lielab::domain::rn>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", &lielab::functions::dexp<lielab::domain::so>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", &lielab::functions::dexp<lielab::domain::sp>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", &lielab::functions::dexp<lielab::domain::su>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &lielab::functions::dexpinv<lielab::domain::gl>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &lielab::functions::dexpinv<lielab::domain::rn>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &lielab::functions::dexpinv<lielab::domain::so>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &lielab::functions::dexpinv<lielab::domain::sp>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &lielab::functions::dexpinv<lielab::domain::su>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("left_exp_default", &lielab::functions::left_exp_default, "Default left exponential action.");

    /*!
    * Begin content for the "kinematics" submodule.
    */
    py::module m_kinematics = m.def_submodule("kinematics", "The kinematics submodule.");
    m_kinematics.def("quaternions", &lielab::kinematics::quaternions, "The quaternion kinematic equation of motion."); // TODO: Change this to use the actual quaternion from domain
    m_kinematics.def("eanglespace123", &lielab::kinematics::eanglespace123, "The eanglespace123 equation of motion.");
    m_kinematics.def("eanglespace231", &lielab::kinematics::eanglespace231, "The eanglespace231 equation of motion.");
    m_kinematics.def("eanglespace312", &lielab::kinematics::eanglespace312, "The eanglespace312 equation of motion.");
    m_kinematics.def("eanglespace132", &lielab::kinematics::eanglespace132, "The eanglespace132 equation of motion.");
    m_kinematics.def("eanglespace213", &lielab::kinematics::eanglespace213, "The eanglespace213 equation of motion.");
    m_kinematics.def("eanglespace321", &lielab::kinematics::eanglespace321, "The eanglespace321 equation of motion.");
    m_kinematics.def("eanglespace121", &lielab::kinematics::eanglespace121, "The eanglespace121 equation of motion.");
    m_kinematics.def("eanglespace131", &lielab::kinematics::eanglespace131, "The eanglespace131 equation of motion.");
    m_kinematics.def("eanglespace212", &lielab::kinematics::eanglespace212, "The eanglespace212 equation of motion.");
    m_kinematics.def("eanglespace232", &lielab::kinematics::eanglespace232, "The eanglespace232 equation of motion.");
    m_kinematics.def("eanglespace313", &lielab::kinematics::eanglespace313, "The eanglespace313 equation of motion.");
    m_kinematics.def("eanglespace323", &lielab::kinematics::eanglespace323, "The eanglespace323 equation of motion.");
    m_kinematics.def("eanglebody123", &lielab::kinematics::eanglebody123, "The eanglebody123 equation of motion.");
    m_kinematics.def("eanglebody231", &lielab::kinematics::eanglebody231, "The eanglebody231 equation of motion.");
    m_kinematics.def("eanglebody312", &lielab::kinematics::eanglebody312, "The eanglebody312 equation of motion.");
    m_kinematics.def("eanglebody132", &lielab::kinematics::eanglebody132, "The eanglebody132 equation of motion.");
    m_kinematics.def("eanglebody213", &lielab::kinematics::eanglebody213, "The eanglebody213 equation of motion.");
    m_kinematics.def("eanglebody321", &lielab::kinematics::eanglebody321, "The eanglebody321 equation of motion.");
    m_kinematics.def("eanglebody121", &lielab::kinematics::eanglebody121, "The eanglebody121 equation of motion.");
    m_kinematics.def("eanglebody131", &lielab::kinematics::eanglebody131, "The eanglebody131 equation of motion.");
    m_kinematics.def("eanglebody212", &lielab::kinematics::eanglebody212, "The eanglebody212 equation of motion.");
    m_kinematics.def("eanglebody232", &lielab::kinematics::eanglebody232, "The eanglebody232 equation of motion.");
    m_kinematics.def("eanglebody313", &lielab::kinematics::eanglebody313, "The eanglebody313 equation of motion.");
    m_kinematics.def("eanglebody323", &lielab::kinematics::eanglebody323, "The eanglebody323 equation of motion.");

    /*!
    * Begin content for the "optim" submodule.
    */
    py::module m_optim = m.def_submodule("optim", "The optim submodule.");

    py::class_<lielab::optim::opt_golden>(m_optim, "opt_golden")
        .def(py::init())
        .def("init", &lielab::optim::opt_golden::init)
        .def("step", &lielab::optim::opt_golden::step)
        .def_readwrite("iterations", &lielab::optim::opt_golden::iterations)
        .def_readwrite("max_iterations", &lielab::optim::opt_golden::max_iterations)
        .def_readwrite("num_objective_evals", &lielab::optim::opt_golden::num_objective_evals)
        .def_readwrite("num_jacobian_evals", &lielab::optim::opt_golden::num_jacobian_evals)
        .def_readwrite("num_hessian_evals", &lielab::optim::opt_golden::num_hessian_evals)
        .def_readwrite("algo_status", &lielab::optim::opt_golden::algo_status)
        .def_readwrite("success", &lielab::optim::opt_golden::success)
        .def_readwrite("tolerance", &lielab::optim::opt_golden::tolerance)
        .def_readwrite("val_objective", &lielab::optim::opt_golden::val_objective)
        .def_readwrite("val_jacobian", &lielab::optim::opt_golden::val_jacobian)
        .def_readwrite("val_hessian", &lielab::optim::opt_golden::val_hessian)
        .def_readwrite("x0", &lielab::optim::opt_golden::x0)
        .def_readwrite("lower", &lielab::optim::opt_golden::lower)
        .def_readwrite("upper", &lielab::optim::opt_golden::upper)
        .def_readwrite("tau", &lielab::optim::opt_golden::tau)
        .def_readwrite("_f1", &lielab::optim::opt_golden::_f1)
        .def_readwrite("_f2", &lielab::optim::opt_golden::_f2)
        .def_readwrite("_X", &lielab::optim::opt_golden::_X)
        .def_readwrite("_X1", &lielab::optim::opt_golden::_X1)
        .def_readwrite("_X2", &lielab::optim::opt_golden::_X2)
        .def_readwrite("_A", &lielab::optim::opt_golden::_A)
        .def_readwrite("_B", &lielab::optim::opt_golden::_B);
    
    py::class_<lielab::optim::hnewton>(m_optim, "hnewton")
        .def(py::init())
        .def("init", &lielab::optim::hnewton::init)
        .def("step0", &lielab::optim::hnewton::step0)
        .def("step1", &lielab::optim::hnewton::step1)
        .def_readwrite("lower", &lielab::optim::hnewton::lower)
        .def_readwrite("upper", &lielab::optim::hnewton::upper)
        .def_readwrite("algo_status", &lielab::optim::hnewton::algo_status)
        .def_readwrite("f", &lielab::optim::hnewton::f)
        .def_readwrite("fnext", &lielab::optim::hnewton::fnext)
        .def_readwrite("m", &lielab::optim::hnewton::m)
        .def_readwrite("mnext", &lielab::optim::hnewton::mnext);
    
    // py::class_<lielab::optim::search_linearx>(m_optim, "search_linearx")
    //     .def(py::init())
    //     .def("init", &lielab::optim::search_linearx::init)
    //     .def("step", &lielab::optim::search_linearx::step)
    //     .def_readwrite("iterations", &lielab::optim::search_linearx::iterations)
    //     .def_readwrite("max_iterations", &lielab::optim::search_linearx::max_iterations)
    //     .def_readwrite("num_objective_evals", &lielab::optim::search_linearx::num_objective_evals)
    //     .def_readwrite("num_jacobian_evals", &lielab::optim::search_linearx::num_jacobian_evals)
    //     .def_readwrite("num_hessian_evals", &lielab::optim::search_linearx::num_hessian_evals)
    //     .def_readwrite("algo_status", &lielab::optim::search_linearx::algo_status)
    //     .def_readwrite("success", &lielab::optim::search_linearx::success)
    //     .def_readwrite("tolerance", &lielab::optim::search_linearx::tolerance)
    //     .def_readwrite("val_objective", &lielab::optim::search_linearx::val_objective)
    //     .def_readwrite("val_jacobian", &lielab::optim::search_linearx::val_jacobian)
    //     .def_readwrite("val_hessian", &lielab::optim::search_linearx::val_hessian)
    //     .def_readwrite("_dx", &lielab::optim::search_linearx::_dx)
    //     .def_readwrite("_x", &lielab::optim::search_linearx::_x)
    //     .def_readwrite("_x1", &lielab::optim::search_linearx::_x1)
    //     .def_readwrite("_x2", &lielab::optim::search_linearx::_x2)
    //     .def_readwrite("_y1", &lielab::optim::search_linearx::_y1)
    //     .def_readwrite("_y2", &lielab::optim::search_linearx::_y2)
    //     .def_readwrite("k", &lielab::optim::search_linearx::k)
    //     .def_readwrite("_lo", &lielab::optim::search_linearx::_lo)
    //     .def_readwrite("_hi", &lielab::optim::search_linearx::_hi)
    //     .def_readwrite("fdx", &lielab::optim::search_linearx::fdx)
    //     .def_readwrite("lower", &lielab::optim::search_linearx::lower)
    //     .def_readwrite("upper", &lielab::optim::search_linearx::upper);
    

    /*!
    * Begin content for the "transform" submodule.
    */

    py::module m_transform = m.def_submodule("transform", "The transform submodule.");
    m_transform.def("dcm_to_quaternion", &lielab::transform::dcm_to_quaternion, "The dcm_to_quaternion function.");
    m_transform.def("dcm_to_eanglebody123", &lielab::transform::dcm_to_eanglebody123<double>, "The dcm_to_eanglebody123 function.");
    m_transform.def("dcm_to_eanglebody231", &lielab::transform::dcm_to_eanglebody231<double>, "The dcm_to_eanglebody231 function.");
    m_transform.def("dcm_to_eanglebody312", &lielab::transform::dcm_to_eanglebody312<double>, "The dcm_to_eanglebody312 function.");
    m_transform.def("dcm_to_eanglebody132", &lielab::transform::dcm_to_eanglebody132<double>, "The dcm_to_eanglebody132 function.");
    m_transform.def("dcm_to_eanglebody213", &lielab::transform::dcm_to_eanglebody213<double>, "The dcm_to_eanglebody213 function.");
    m_transform.def("dcm_to_eanglebody321", &lielab::transform::dcm_to_eanglebody321<double>, "The dcm_to_eanglebody321 function.");
    m_transform.def("dcm_to_eanglebody121", &lielab::transform::dcm_to_eanglebody121<double>, "The dcm_to_eanglebody121 function.");
    m_transform.def("dcm_to_eanglebody131", &lielab::transform::dcm_to_eanglebody131<double>, "The dcm_to_eanglebody131 function.");
    m_transform.def("dcm_to_eanglebody212", &lielab::transform::dcm_to_eanglebody212<double>, "The dcm_to_eanglebody212 function.");
    m_transform.def("dcm_to_eanglebody232", &lielab::transform::dcm_to_eanglebody232<double>, "The dcm_to_eanglebody232 function.");
    m_transform.def("dcm_to_eanglebody313", &lielab::transform::dcm_to_eanglebody313<double>, "The dcm_to_eanglebody313 function.");
    m_transform.def("dcm_to_eanglebody323", &lielab::transform::dcm_to_eanglebody323<double>, "The dcm_to_eanglebody323 function.");
    m_transform.def("dcm_to_eanglespace123", &lielab::transform::dcm_to_eanglespace123<double>, "The dcm_to_eanglespace123 function.");
    m_transform.def("dcm_to_eanglespace231", &lielab::transform::dcm_to_eanglespace231<double>, "The dcm_to_eanglespace231 function.");
    m_transform.def("dcm_to_eanglespace312", &lielab::transform::dcm_to_eanglespace312<double>, "The dcm_to_eanglespace312 function.");
    m_transform.def("dcm_to_eanglespace132", &lielab::transform::dcm_to_eanglespace132<double>, "The dcm_to_eanglespace132 function.");
    m_transform.def("dcm_to_eanglespace213", &lielab::transform::dcm_to_eanglespace213<double>, "The dcm_to_eanglespace213 function.");
    m_transform.def("dcm_to_eanglespace321", &lielab::transform::dcm_to_eanglespace321<double>, "The dcm_to_eanglespace321 function.");
    m_transform.def("dcm_to_eanglespace121", &lielab::transform::dcm_to_eanglespace121<double>, "The dcm_to_eanglespace121 function.");
    m_transform.def("dcm_to_eanglespace131", &lielab::transform::dcm_to_eanglespace131<double>, "The dcm_to_eanglespace131 function.");
    m_transform.def("dcm_to_eanglespace212", &lielab::transform::dcm_to_eanglespace212<double>, "The dcm_to_eanglespace212 function.");
    m_transform.def("dcm_to_eanglespace232", &lielab::transform::dcm_to_eanglespace232<double>, "The dcm_to_eanglespace232 function.");
    m_transform.def("dcm_to_eanglespace313", &lielab::transform::dcm_to_eanglespace313<double>, "The dcm_to_eanglespace313 function.");
    m_transform.def("dcm_to_eanglespace323", &lielab::transform::dcm_to_eanglespace323<double>, "The dcm_to_eanglespace323 function.");
    m_transform.def("eanglebody123_to_dcm", &lielab::transform::eanglebody123_to_dcm<double>, "The eanglebody123_to_dcm function.");
    m_transform.def("eanglebody231_to_dcm", &lielab::transform::eanglebody231_to_dcm<double>, "The eanglebody231_to_dcm function.");
    m_transform.def("eanglebody312_to_dcm", &lielab::transform::eanglebody312_to_dcm<double>, "The eanglebody312_to_dcm function.");
    m_transform.def("eanglebody132_to_dcm", &lielab::transform::eanglebody132_to_dcm<double>, "The eanglebody132_to_dcm function.");
    m_transform.def("eanglebody213_to_dcm", &lielab::transform::eanglebody213_to_dcm<double>, "The eanglebody213_to_dcm function.");
    m_transform.def("eanglebody321_to_dcm", &lielab::transform::eanglebody321_to_dcm<double>, "The eanglebody321_to_dcm function.");
    m_transform.def("eanglebody121_to_dcm", &lielab::transform::eanglebody121_to_dcm<double>, "The eanglebody121_to_dcm function.");
    m_transform.def("eanglebody131_to_dcm", &lielab::transform::eanglebody131_to_dcm<double>, "The eanglebody131_to_dcm function.");
    m_transform.def("eanglebody212_to_dcm", &lielab::transform::eanglebody212_to_dcm<double>, "The eanglebody212_to_dcm function.");
    m_transform.def("eanglebody232_to_dcm", &lielab::transform::eanglebody232_to_dcm<double>, "The eanglebody232_to_dcm function.");
    m_transform.def("eanglebody313_to_dcm", &lielab::transform::eanglebody313_to_dcm<double>, "The eanglebody313_to_dcm function.");
    m_transform.def("eanglebody323_to_dcm", &lielab::transform::eanglebody323_to_dcm<double>, "The eanglebody323_to_dcm function.");
    m_transform.def("eanglespace123_to_dcm", &lielab::transform::eanglespace123_to_dcm<double>, "The eanglespace123_to_dcm function.");
    m_transform.def("eanglespace231_to_dcm", &lielab::transform::eanglespace231_to_dcm<double>, "The eanglespace231_to_dcm function.");
    m_transform.def("eanglespace312_to_dcm", &lielab::transform::eanglespace312_to_dcm<double>, "The eanglespace312_to_dcm function.");
    m_transform.def("eanglespace132_to_dcm", &lielab::transform::eanglespace132_to_dcm<double>, "The eanglespace132_to_dcm function.");
    m_transform.def("eanglespace213_to_dcm", &lielab::transform::eanglespace213_to_dcm<double>, "The eanglespace213_to_dcm function.");
    m_transform.def("eanglespace321_to_dcm", &lielab::transform::eanglespace321_to_dcm<double>, "The eanglespace321_to_dcm function.");
    m_transform.def("eanglespace121_to_dcm", &lielab::transform::eanglespace121_to_dcm<double>, "The eanglespace121_to_dcm function.");
    m_transform.def("eanglespace131_to_dcm", &lielab::transform::eanglespace131_to_dcm<double>, "The eanglespace131_to_dcm function.");
    m_transform.def("eanglespace212_to_dcm", &lielab::transform::eanglespace212_to_dcm<double>, "The eanglespace212_to_dcm function.");
    m_transform.def("eanglespace232_to_dcm", &lielab::transform::eanglespace232_to_dcm<double>, "The eanglespace232_to_dcm function.");
    m_transform.def("eanglespace313_to_dcm", &lielab::transform::eanglespace313_to_dcm<double>, "The eanglespace313_to_dcm function.");
    m_transform.def("eanglespace323_to_dcm", &lielab::transform::eanglespace323_to_dcm<double>, "The eanglespace323_to_dcm function.");
    m_transform.def("quaternion_to_dcm", &lielab::transform::quaternion_to_dcm, "The quaternion_to_dcm function.");

    py::module m_topos = m.def_submodule("topos", "The topos submodule.");

    m_topos.def("Ad", &lielab::topos::Ad);
    m_topos.def("dexp", &lielab::topos::dexp, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_topos.def("dexpinv", &lielab::topos::dexpinv);
    m_topos.def("exp", &lielab::topos::exp);
    m_topos.def("log", &lielab::topos::log);
    
    // py::enum_<lielab::topos::RKTYPE>(m_topos, "RKTYPE")
    //     .value("RKTYPE_NONE", lielab::topos::RKTYPE::RKTYPE_NONE)
    //     .value("RKTYPE_EXPLICIT", lielab::topos::RKTYPE::RKTYPE_EXPLICIT)
    //     .value("RKTYPE_IMPLICIT", lielab::topos::RKTYPE::RKTYPE_IMPLICIT);
    
    py::class_<lielab::topos::IntegralCurve>(m_topos, "IntegralCurve")
        .def(py::init<>())
        .def(py::init<lielab::topos::IntegralCurve>())
        .def(py::init<size_t>())
        .def_readwrite("chunk", &lielab::topos::IntegralCurve::chunk)
        .def_readwrite("length", &lielab::topos::IntegralCurve::length)
        .def_readwrite("num_eoms", &lielab::topos::IntegralCurve::num_eoms)
        .def_readwrite("chs", &lielab::topos::IntegralCurve::chs)
        .def_readwrite("t", &lielab::topos::IntegralCurve::t)
        .def_readwrite("y", &lielab::topos::IntegralCurve::y);
    
    py::class_<lielab::topos::TSOutput>(m_topos, "TSOutput")
        .def(py::init())
        .def(py::init<lielab::domain::hmlie>())
        .def(py::init<lielab::domain::hmlie, lielab::domain::hmlie, double>())
        .def_readonly("low", &lielab::topos::TSOutput::low)
        .def_readonly("high", &lielab::topos::TSOutput::high)
        .def_readonly("error", &lielab::topos::TSOutput::error);

    py::enum_<lielab::topos::RKMETHOD>(m_topos, "RKMETHOD")
        .value("E1", lielab::topos::RKMETHOD::E1)
        .value("RK45", lielab::topos::RKMETHOD::RK45);
    
    py::class_<lielab::topos::TimeStepper>(m_topos, "TimeStepper")
        .def_readwrite("algo_status", &lielab::topos::TimeStepper::algo_status)
        .def_readwrite("iterations", &lielab::topos::TimeStepper::iterations)
        .def_readwrite("max_iterations", &lielab::topos::TimeStepper::max_iterations)
        .def_readwrite("success", &lielab::topos::TimeStepper::success)
        .def_readonly("A", &lielab::topos::TimeStepper::A)
        .def_readonly("B", &lielab::topos::TimeStepper::B)
        .def_readonly("Bhat", &lielab::topos::TimeStepper::Bhat)
        .def_readonly("C", &lielab::topos::TimeStepper::C)
        .def_readonly("n", &lielab::topos::TimeStepper::n)
        .def_readonly("order", &lielab::topos::TimeStepper::order)
        .def_readonly("variable_step", &lielab::topos::TimeStepper::variable_step)
        .def(py::init())
        .def(py::init<lielab::topos::RKMETHOD>())
        .def("init", &lielab::topos::TimeStepper::init)
        .def("step", &lielab::topos::TimeStepper::step);

    py::class_<lielab::topos::MuntheKaas>(m_topos, "MuntheKaas")
        .def_readwrite("algo_status", &lielab::topos::MuntheKaas::algo_status)
        .def_readwrite("iterations", &lielab::topos::TimeStepper::iterations)
        .def_readwrite("max_iterations", &lielab::topos::TimeStepper::max_iterations)
        .def_readwrite("success", &lielab::topos::TimeStepper::success)
        .def_readonly("_KK", &lielab::topos::MuntheKaas::_KK)
        .def_readonly("_U", &lielab::topos::MuntheKaas::_U)
        .def_readonly("A", &lielab::topos::MuntheKaas::A)
        .def_readonly("B", &lielab::topos::MuntheKaas::B)
        .def_readonly("Bhat", &lielab::topos::MuntheKaas::Bhat)
        .def_readonly("C", &lielab::topos::MuntheKaas::C)
        .def_readonly("n", &lielab::topos::MuntheKaas::n)
        .def_readonly("order", &lielab::topos::MuntheKaas::order)
        .def_readonly("variable_step", &lielab::topos::MuntheKaas::variable_step)
        .def(py::init())
        .def(py::init<lielab::topos::RKMETHOD>())
        .def("init", &lielab::topos::MuntheKaas::init)
        .def("step_0", &lielab::topos::MuntheKaas::step_0)
        .def("step_1", &lielab::topos::MuntheKaas::step_1)
        .def("set_dy", &lielab::topos::MuntheKaas::set_dy)
        .def("postprocess", &lielab::topos::MuntheKaas::postprocess)
        .def_readwrite("next_t", &lielab::topos::MuntheKaas::next_t)
        .def_readwrite("next_y", &lielab::topos::MuntheKaas::next_y)
        .def("__repr__",
        [](const lielab::topos::MuntheKaas & self)
        {
            return "<lielab.topos.MuntheKaas>";
        })
        .def("__str__",
        [](const lielab::topos::MuntheKaas & self)
        {
            return "<lielab.topos.MuntheKaas>";
        });
    
    py::class_<lielab::topos::Flow>(m_topos, "Flow")
        .def(py::init())
        .def("init", &lielab::topos::Flow::init)
        .def("step0", &lielab::topos::Flow::step0)
        .def("step", &lielab::topos::Flow::step)
        .def("stepE", &lielab::topos::Flow::stepE)
        .def("postprocess", &lielab::topos::Flow::postprocess)
        .def("new_step_size", &lielab::topos::Flow::new_step_size)
        .def_readwrite("algo_status", &lielab::topos::Flow::algo_status)
        .def_readwrite("_ynext", &lielab::topos::Flow::_ynext)
        .def_readwrite("iterations", &lielab::topos::Flow::iterations)
        .def_readwrite("small", &lielab::topos::Flow::small)
        .def_readwrite("large", &lielab::topos::Flow::large)
        .def_readwrite("pessimist", &lielab::topos::Flow::pessimist)
        .def_readwrite("accept", &lielab::topos::Flow::accept)
        .def_readwrite("tol", &lielab::topos::Flow::tol)
        .def_readwrite("dt_min", &lielab::topos::Flow::dt_min)
        .def_readwrite("dt_max", &lielab::topos::Flow::dt_max)
        .def_readwrite("dt", &lielab::topos::Flow::dt)
        .def_readwrite("default_local", &lielab::topos::Flow::default_local)
        .def_readwrite("default_global", &lielab::topos::Flow::default_global)
        .def_readwrite("variable_time_step", &lielab::topos::Flow::variable_time_step)
    //     .def_readwrite("new_exact", &lielab::topos::Flow::new_exact)
    //     .def_readwrite("num_step", &lielab::topos::Flow::num_step)
        .def_readwrite("_out", &lielab::topos::Flow::_out)
        .def_readwrite("_dt", &lielab::topos::Flow::_dt)
        .def_readwrite("stepper", &lielab::topos::Flow::stepper)
        .def("__repr__",
        [](const lielab::topos::Flow & self)
        {
            return "<lielab.topos.Flow>";
        })
        .def("__str__",
        [](const lielab::topos::Flow & self)
        {
            return "<lielab.topos.Flow>";
        });


    // /*!
    // * Begin content for the "transform" submodule.
    // */
    // py::module m_transform = m.def_submodule("transform", "The transform submodule.");
    // m_transform.def("dcm_to_quaternion", &lielab::transform::dcm_to_quaternion, "The dcm to quaternion transform function.");
    // m_transform.def("eanglebody123_to_dcm", &lielab::transform::eanglebody123_to_dcm, "Euler angle body-123 rotation.");
    // m_transform.def("eanglebody231_to_dcm", &lielab::transform::eanglebody231_to_dcm, "Euler angle body-231 rotation.");
    // m_transform.def("eanglebody312_to_dcm", &lielab::transform::eanglebody312_to_dcm, "Euler angle body-312 rotation.");
    // m_transform.def("eanglebody132_to_dcm", &lielab::transform::eanglebody132_to_dcm, "Euler angle body-132 rotation.");
    // m_transform.def("eanglebody213_to_dcm", &lielab::transform::eanglebody213_to_dcm, "Euler angle body-213 rotation.");
    // m_transform.def("eanglebody321_to_dcm", &lielab::transform::eanglebody321_to_dcm, "Euler angle body-321 rotation.");
    // m_transform.def("eanglebody121_to_dcm", &lielab::transform::eanglebody121_to_dcm, "Euler angle body-121 rotation.");
    // m_transform.def("eanglebody131_to_dcm", &lielab::transform::eanglebody131_to_dcm, "Euler angle body-131 rotation.");
    // m_transform.def("eanglebody212_to_dcm", &lielab::transform::eanglebody212_to_dcm, "Euler angle body-212 rotation.");
    // m_transform.def("eanglebody232_to_dcm", &lielab::transform::eanglebody232_to_dcm, "Euler angle body-232 rotation.");
    // m_transform.def("eanglebody313_to_dcm", &lielab::transform::eanglebody313_to_dcm, "Euler angle body-313 rotation.");
    // m_transform.def("eanglebody323_to_dcm", &lielab::transform::eanglebody323_to_dcm, "Euler angle body-323 rotation.");
    // m_transform.def("eanglespace123_to_dcm", &lielab::transform::eanglespace123_to_dcm, "Euler angle space-123 rotation.");
    // m_transform.def("eanglespace231_to_dcm", &lielab::transform::eanglespace231_to_dcm, "Euler angle space-231 rotation.");
    // m_transform.def("eanglespace312_to_dcm", &lielab::transform::eanglespace312_to_dcm, "Euler angle space-312 rotation.");
    // m_transform.def("eanglespace132_to_dcm", &lielab::transform::eanglespace132_to_dcm, "Euler angle space-132 rotation.");
    // m_transform.def("eanglespace213_to_dcm", &lielab::transform::eanglespace213_to_dcm, "Euler angle space-213 rotation.");
    // m_transform.def("eanglespace321_to_dcm", &lielab::transform::eanglespace321_to_dcm, "Euler angle space-321 rotation.");
    // m_transform.def("eanglespace121_to_dcm", &lielab::transform::eanglespace121_to_dcm, "Euler angle space-121 rotation.");
    // m_transform.def("eanglespace131_to_dcm", &lielab::transform::eanglespace131_to_dcm, "Euler angle space-131 rotation.");
    // m_transform.def("eanglespace212_to_dcm", &lielab::transform::eanglespace212_to_dcm, "Euler angle space-212 rotation.");
    // m_transform.def("eanglespace232_to_dcm", &lielab::transform::eanglespace232_to_dcm, "Euler angle space-232 rotation.");
    // m_transform.def("eanglespace313_to_dcm", &lielab::transform::eanglespace313_to_dcm, "Euler angle space-313 rotation.");
    // m_transform.def("eanglespace323_to_dcm", &lielab::transform::eanglespace323_to_dcm, "Euler angle space-323 rotation.");
    // m_transform.def("gibbs_to_dcm", &lielab::transform::gibbs_to_dcm, "Gibbs to DCM.");
    // m_transform.def("dcm_to_gibbs", &lielab::transform::dcm_to_gibbs, "DCM to Gibbs.");
    // m_transform.def("quaternion_to_dcm", &lielab::transform::quaternion_to_dcm, "The quaternion to dcm transform function.");
}
