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
    * Bindings for Lielab::domain::gl
    */

    py::class_<Lielab::domain::gl>(m_domain, "gl")
        .def_readwrite("_data", &Lielab::domain::gl::_data)
        .def_readonly_static("abelian", &Lielab::domain::gl::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXd &>())
        .def("basis", &Lielab::domain::gl::basis)
        .def("project", &Lielab::domain::gl::project)
        .def("get_dimension", &Lielab::domain::gl::get_dimension)
        .def("get_vector", &Lielab::domain::gl::get_vector)
        .def("get_ados_representation", &Lielab::domain::gl::get_ados_representation)
        .def("set_vector", &Lielab::domain::gl::set_vector)
        .def("__call__", [](Lielab::domain::gl & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const Lielab::domain::gl & self, const size_t index1, const size_t index2)
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
        .def("__repr__", [](const Lielab::domain::gl & self)
        {
            return "<lielab.domain.gl>";
        })
        .def("__str__", [](const Lielab::domain::gl & self)
        {
            return matstr(self._data);
        });

    /*!
    * Bindings for Lielab::domain::rn
    */

    py::class_<Lielab::domain::rn>(m_domain, "rn")
        .def_readwrite("_data", &Lielab::domain::rn::_data)
        .def_readwrite("shape", &Lielab::domain::rn::shape)
        .def_readonly_static("abelian", &Lielab::domain::rn::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXd &>())
        .def("basis", &Lielab::domain::rn::basis)
        .def("project", &Lielab::domain::rn::project)
        .def("get_dimension", &Lielab::domain::rn::get_dimension)
        .def("get_vector", &Lielab::domain::rn::get_vector)
        .def("get_ados_representation", &Lielab::domain::rn::get_ados_representation)
        .def("set_vector", &Lielab::domain::rn::set_vector)
        .def("__call__", [](const Lielab::domain::rn & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const Lielab::domain::rn & self, const size_t index1, const size_t index2)
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
        .def("__repr__", [](const Lielab::domain::rn & self)
        {
            return "<lielab.domain.rn>";
        })
        .def("__str__", [](const Lielab::domain::rn & self)
        {
            return matstr(self._data);
        });
    
    /*!
    * Bindings for Lielab::domain::se
    */

    py::class_<Lielab::domain::se>(m_domain, "se")
        .def_readwrite("_data", &Lielab::domain::se::_data)
        .def_readwrite("shape", &Lielab::domain::se::shape)
        .def_readonly_static("abelian", &Lielab::domain::se::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXd &>())
        .def("basis", &Lielab::domain::se::basis)
        // .def("project", &Lielab::domain::se::project)
        .def("get_dimension", &Lielab::domain::se::get_dimension)
        .def("get_vector", &Lielab::domain::se::get_vector)
        .def("get_ados_representation", &Lielab::domain::se::get_ados_representation)
        .def("set_vector", &Lielab::domain::se::set_vector)
        .def("__call__", [](const Lielab::domain::se & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const Lielab::domain::se & self, const size_t index1, const size_t index2)
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
        .def("__repr__", [](const Lielab::domain::se & self)
        {
            return "<lielab.domain.se>";
        })
        .def("__str__", [](const Lielab::domain::se & self)
        {
            return matstr(self._data);
        });


    /*!
    * Bindings for Lielab::domain::so
    */

    py::class_<Lielab::domain::so>(m_domain, "so")
        .def_readwrite("_data", &Lielab::domain::so::_data)
        .def_readwrite("shape", &Lielab::domain::so::shape)
        .def_readonly_static("abelian", &Lielab::domain::so::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXd &>())
        .def("basis", &Lielab::domain::so::basis)
        .def("project", &Lielab::domain::so::project)
        .def("get_dimension", &Lielab::domain::so::get_dimension)
        .def("get_vector", &Lielab::domain::so::get_vector)
        .def("get_ados_representation", &Lielab::domain::so::get_ados_representation)
        .def("set_vector", &Lielab::domain::so::set_vector)
        .def("__call__", [](const Lielab::domain::so & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const Lielab::domain::so & self, const size_t index1, const size_t index2)
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
        .def("__repr__", [](const Lielab::domain::so & self)
        {
            return "<lielab.domain.so>";
        })
        .def("__str__", [](const Lielab::domain::so & self)
        {
            return matstr(self._data);
        });

    
    /*!
    * Bindings for Lielab::domain::sp
    */
    
    py::class_<Lielab::domain::sp>(m_domain, "sp")
        .def_readwrite("_data", &Lielab::domain::sp::_data)
        .def_readwrite("shape", &Lielab::domain::sp::shape)
        .def_readonly_static("abelian", &Lielab::domain::sp::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXd &>())
        .def("basis", &Lielab::domain::sp::basis)
        .def("project", &Lielab::domain::sp::project)
        .def("get_dimension", &Lielab::domain::sp::get_dimension)
        .def("get_vector", &Lielab::domain::sp::get_vector)
        .def("get_ados_representation", &Lielab::domain::sp::get_ados_representation)
        .def("set_vector", &Lielab::domain::sp::set_vector)
        .def("__call__", [](const Lielab::domain::sp & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const Lielab::domain::sp & self, const size_t index1, const size_t index2)
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
        .def("__repr__", [](const Lielab::domain::sp & self)
        {
            return "<lielab.domain.sp>";
        })
        .def("__str__", [](const Lielab::domain::sp & self)
        {
            return matstr(self._data);
        });
    
    
    /*!
    * Bindings for Lielab::domain::su
    */
    
    py::class_<Lielab::domain::su>(m_domain, "su")
        .def_readwrite("_data", &Lielab::domain::su::_data)
        .def_readwrite("shape", &Lielab::domain::su::shape)
        .def_readonly_static("abelian", &Lielab::domain::su::abelian)
        .def(py::init<>())
        .def(py::init<const size_t>())
        .def(py::init<const Eigen::MatrixXcd &>())
        .def("basis", &Lielab::domain::su::basis)
        // .def("project", &Lielab::domain::su::project) // TODO: Implement project
        .def("get_dimension", &Lielab::domain::su::get_dimension)
        .def("get_vector", &Lielab::domain::su::get_vector)
        .def("get_ados_representation", &Lielab::domain::su::get_ados_representation)
        .def("set_vector", &Lielab::domain::su::set_vector)
        .def("__call__", [](const Lielab::domain::su & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const Lielab::domain::su & self, const size_t index1, const size_t index2)
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
        .def("__repr__", [](const Lielab::domain::su & self)
        {
            return "<lielab.domain.su>";
        })
        .def("__str__", [](const Lielab::domain::su & self)
        {
            return matstr(self._data);
        });
    
    /*!
     * Lie Groups
     */
    
    py::class_<Lielab::domain::GL>(m_domain, "GL")
        .def_readwrite("_data", &Lielab::domain::GL::_data)
        .def_readonly_static("abelian", &Lielab::domain::GL::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXd>())
        .def("project", &Lielab::domain::GL::project)
        .def("get_dimension", &Lielab::domain::GL::get_dimension)
        .def("get_ados_representation", &Lielab::domain::GL::get_ados_representation)
        .def("inverse", &Lielab::domain::GL::inverse)
        .def("serialize", &Lielab::domain::GL::serialize)
        .def("unserialize", &Lielab::domain::GL::unserialize)
        .def("__call__", [](const Lielab::domain::GL & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const Lielab::domain::GL & self)
        {
            return "<lielab.domain.GL>";
        })
        .def("__str__", [](const Lielab::domain::GL & self)
        {
            return matstr(self._data);
        });

    py::class_<Lielab::domain::RN>(m_domain, "RN")
        .def_readwrite("_data", &Lielab::domain::RN::_data)
        .def_readwrite("shape", &Lielab::domain::RN::shape)
        .def_readonly_static("abelian", &Lielab::domain::RN::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXd>())
        .def("project", &Lielab::domain::RN::project)
        .def("get_dimension", &Lielab::domain::RN::get_dimension)
        .def("get_ados_representation", &Lielab::domain::RN::get_ados_representation)
        .def("inverse", &Lielab::domain::RN::inverse)
        .def("serialize", &Lielab::domain::RN::serialize)
        .def("unserialize", &Lielab::domain::RN::unserialize)
        .def("__call__", [](const Lielab::domain::RN & self, const size_t index)
        {
            return self(index);
        })
        .def("__call__", [](const Lielab::domain::RN & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const Lielab::domain::RN & self)
        {
            return "<lielab.domain.RN>";
        })
        .def("__str__", [](const Lielab::domain::RN & self)
        {
            return matstr(self._data);
        });

    py::class_<Lielab::domain::SE>(m_domain, "SE")
        .def_readwrite("_data", &Lielab::domain::SE::_data)
        .def_readwrite("shape", &Lielab::domain::SE::shape)
        .def_readonly_static("abelian", &Lielab::domain::SE::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXd>())
        // .def("project", &Lielab::domain::SE::project)
        .def("get_dimension", &Lielab::domain::SE::get_dimension)
        .def("get_ados_representation", &Lielab::domain::SE::get_ados_representation)
        .def("inverse", &Lielab::domain::SE::inverse)
        .def("serialize", &Lielab::domain::SE::serialize)
        .def("unserialize", &Lielab::domain::SE::unserialize)
        .def("__call__", [](const Lielab::domain::SE & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const Lielab::domain::SE & self)
        {
            return "<lielab.domain.SE>";
        })
        .def("__str__", [](const Lielab::domain::SE & self)
        {
            return matstr(self._data);
        });

    py::class_<Lielab::domain::SO>(m_domain, "SO")
        .def_readwrite("_data", &Lielab::domain::SO::_data)
        .def_readwrite("shape", &Lielab::domain::SO::shape)
        .def_readonly_static("abelian", &Lielab::domain::SO::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXd>())
        .def("project", &Lielab::domain::SO::project)
        .def("get_dimension", &Lielab::domain::SO::get_dimension)
        .def("get_ados_representation", &Lielab::domain::SO::get_ados_representation)
        .def("inverse", &Lielab::domain::SO::inverse)
        .def("serialize", &Lielab::domain::SO::serialize)
        .def("unserialize", &Lielab::domain::SO::unserialize)
        .def("__call__", [](const Lielab::domain::SO & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const Lielab::domain::SO & self)
        {
            return "<lielab.domain.SO>";
        })
        .def("__str__", [](const Lielab::domain::SO & self)
        {
            return matstr(self._data);
        });
    

    py::class_<Lielab::domain::SP>(m_domain, "SP")
        .def_readwrite("_data", &Lielab::domain::SP::_data)
        .def_readwrite("shape", &Lielab::domain::SP::shape)
        .def_readonly_static("abelian", &Lielab::domain::SP::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXd>())
        // .def("project", &Lielab::domain::SP::project) // TODO:
        .def("get_dimension", &Lielab::domain::SP::get_dimension)
        .def("get_ados_representation", &Lielab::domain::SP::get_ados_representation)
        .def("inverse", &Lielab::domain::SP::inverse)
        .def("serialize", &Lielab::domain::SP::serialize)
        .def("unserialize", &Lielab::domain::SP::unserialize)
        .def("__call__", [](const Lielab::domain::SP & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const Lielab::domain::SP & self)
        {
            return "<lielab.domain.SP>";
        })
        .def("__str__", [](const Lielab::domain::SP & self)
        {
            return matstr(self._data);
        });

    py::class_<Lielab::domain::SU>(m_domain, "SU")
        .def_readwrite("_data", &Lielab::domain::SU::_data)
        .def_readwrite("shape", &Lielab::domain::SU::shape)
        .def_readonly_static("abelian", &Lielab::domain::SU::abelian)
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<Eigen::MatrixXcd>())
        .def_static("Quaternion", py::overload_cast<>(&Lielab::domain::SU::Quaternion))
        .def_static("Quaternion", py::overload_cast<const double, const double, const double, const double>(&Lielab::domain::SU::Quaternion))
        // .def("project", &Lielab::domain::SU::project) // TODO:
        .def("get_dimension", &Lielab::domain::SU::get_dimension)
        .def("get_ados_representation", &Lielab::domain::SU::get_ados_representation)
        .def("inverse", &Lielab::domain::SU::inverse)
        .def("serialize", &Lielab::domain::SU::serialize)
        .def("unserialize", &Lielab::domain::SU::unserialize)
        .def("__call__", [](const Lielab::domain::SU & self, const size_t index1, const size_t index2)
        {
            return self(index1, index2);
        })
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("__repr__", [](const Lielab::domain::SU & self)
        {
            return "<lielab.domain.SU>";
        })
        .def("__str__", [](const Lielab::domain::SU & self)
        {
            return matstr(self._data);
        });

    py::class_<Lielab::domain::CompositeAlgebra>(m_domain, "CompositeAlgebra")
        .def(py::init())
        .def(py::init<std::vector<Lielab::domain::CompositeAlgebra::TYPES>>())
        .def_readwrite("space", &Lielab::domain::CompositeAlgebra::space)
        .def("get_dimension", &Lielab::domain::CompositeAlgebra::get_dimension)
        .def("get_shape", &Lielab::domain::CompositeAlgebra::get_shape)
        .def("get_vector", &Lielab::domain::CompositeAlgebra::get_vector)
        .def("set_vector", &Lielab::domain::CompositeAlgebra::set_vector)
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
        [](const Lielab::domain::CompositeAlgebra & self)
        {
            return "<lielab.domain.CompositeAlgebra>";
        })
        .def("__str__",
        [](const Lielab::domain::CompositeAlgebra & self)
        {
            return self.to_string();
        });
    

    py::class_<Lielab::domain::CompositeManifold>(m_domain, "CompositeManifold")
        .def(py::init())
        .def(py::init<std::vector<Lielab::domain::CompositeManifold::TYPES>>())
        .def_readwrite("space", &Lielab::domain::CompositeManifold::space)
        .def("get_shape", &Lielab::domain::CompositeManifold::get_shape)
        .def("serialize", &Lielab::domain::CompositeManifold::serialize)
        .def("unserialize", &Lielab::domain::CompositeManifold::unserialize)
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def("inverse", &Lielab::domain::CompositeManifold::inverse)
        .def("__repr__",
        [](const Lielab::domain::CompositeManifold & self)
        {
            return "<lielab.domain.CompositeManifold>";
        })
        .def("__str__",
        [](const Lielab::domain::CompositeManifold & self)
        {
            return self.to_string();
        });


    /*!
    * Begin content for the "functions" submodule.
    */

    py::module m_functions = m.def_submodule("functions", "The functions submodule.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::gl>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::rn>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::se>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::so>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::sp>, "The pair function.");
    m_functions.def("pair", &Lielab::functions::pair<Lielab::domain::su>, "The pair function.");
    // m_functions.def("copair", &Lielab::functions::copair<Lielab::domain::gl>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<Lielab::domain::rn>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<Lielab::domain::so>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<Lielab::domain::sp>, "The copair function.");
    // m_functions.def("copair", &Lielab::functions::copair<Lielab::domain::su>, "The copair function.");
    m_functions.def("factorial", &Lielab::functions::factorial, "The factorial function");
    m_functions.def("Ad", &Lielab::functions::Ad<Lielab::domain::gl>, "The Ad function.");
    m_functions.def("Ad", &Lielab::functions::Ad<Lielab::domain::rn>, "The Ad function.");
    m_functions.def("Ad", &Lielab::functions::Ad<Lielab::domain::se>, "The Ad function.");
    m_functions.def("Ad", &Lielab::functions::Ad<Lielab::domain::so>, "The Ad function.");
    m_functions.def("Ad", &Lielab::functions::Ad<Lielab::domain::sp>, "The Ad function.");
    m_functions.def("Ad", &Lielab::functions::Ad<Lielab::domain::su>, "The Ad function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::gl>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::rn>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::se>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::so>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::sp>, "The commutator function.");
    m_functions.def("commutator", &Lielab::functions::commutator<Lielab::domain::su>, "The commutator function.");
    // m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::gl>, "The cayley1 function.");
    m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::rn>, "The cayley1 function.");
    // m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::se>, "The cayley1 function.");
    m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::so>, "The cayley1 function.");
    m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::sp>, "The cayley1 function.");
    m_functions.def("cayley1", &Lielab::functions::cayley1<Lielab::domain::su>, "The cayley1 function.");
    // m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::gl>, "The cayley2 function.");
    m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::rn>, "The cayley2 function.");
    // m_functions.def("cayley2", &Lielab::functions::cayley1<Lielab::domain::se>, "The cayley2 function.");
    m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::so>, "The cayley2 function.");
    m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::sp>, "The cayley2 function.");
    // m_functions.def("cayley2", &Lielab::functions::cayley2<Lielab::domain::su>, "The cayley2 function.");
    // m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::gl>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::rn>, "The Killing function.");
    // m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::se>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::so>, "The Killing function.");
    m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::sp>, "The Killing function.");
    // m_functions.def("Killing", &Lielab::functions::Killing<Lielab::domain::su>, "The Killing function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::gl>, "The Killingform function.");
    // m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::se>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::rn>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::so>, "The Killingform function.");
    m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::sp>, "The Killingform function.");
    // m_functions.def("Killingform", &Lielab::functions::Killingform<Lielab::domain::su>, "The Killingform function.");
    // m_functions.def("ad", &Lielab::functions::ad<Lielab::domain::gl>, "The ad function.");
    m_functions.def("ad", &Lielab::functions::ad<Lielab::domain::rn>, "The ad function.");
    // m_functions.def("ad", &Lielab::functions::ad<Lielab::domain::se>, "The ad function.");
    m_functions.def("ad", &Lielab::functions::ad<Lielab::domain::so>, "The ad function.");
    m_functions.def("ad", &Lielab::functions::ad<Lielab::domain::sp>, "The ad function.");
    m_functions.def("ad", &Lielab::functions::ad<Lielab::domain::su>, "The ad function.");
    // m_functions.def("coAd", &Lielab::functions::coAd<Lielab::domain::GL>, "The coAd function.");
    // m_functions.def("coAd", &Lielab::functions::coAd<Lielab::domain::RN>, "The coAd function.");
    // m_functions.def("coAd", &Lielab::functions::coAd<Lielab::domain::SO>, "The coAd function.");
    // m_functions.def("coAd", &Lielab::functions::coAd<Lielab::domain::SP>, "The coAd function.");
    // m_functions.def("coAd", &Lielab::functions::coAd<Lielab::domain::SU>, "The coAd function."); // SU needs get_dimension()
    // m_functions.def("coad", &Lielab::functions::coad<Lielab::domain::gl>, "The coad function.");
    // m_functions.def("coad", &Lielab::functions::coad<Lielab::domain::rn>, "The coad function.");
    // m_functions.def("coad", &Lielab::functions::coad<Lielab::domain::so>, "The coad function.");
    // m_functions.def("coad", &Lielab::functions::coad<Lielab::domain::sp>, "The coad function.");
    // m_functions.def("coad", &Lielab::functions::coad<Lielab::domain::su>, "The coad function."); // TODO: su needs basis()
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::gl>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::rn>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::se>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::so>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::sp>, "The exponential function.");
    m_functions.def("exp", &Lielab::functions::exp<Lielab::domain::su>, "The exponential function.");
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::GL>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::RN>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::SE>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::SO>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::SP>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("log", &Lielab::functions::log<Lielab::domain::SU>, "The logarithm function.", py::arg("G"), py::arg("optimize") = false);
    m_functions.def("bernoulli", &Lielab::functions::bernoulli, "The bernoulli function.");
    m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::gl>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::rn>, "The dcayley1inv function.");
    // m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::se>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::so>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::sp>, "The dcayley1inv function.");
    m_functions.def("dcayley1inv", &Lielab::functions::dcayley1inv<Lielab::domain::su>, "The dcayley1inv function.");
    m_functions.def("dexp", &Lielab::functions::dexp<Lielab::domain::gl>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", &Lielab::functions::dexp<Lielab::domain::rn>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", &Lielab::functions::dexp<Lielab::domain::se>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", &Lielab::functions::dexp<Lielab::domain::so>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", &Lielab::functions::dexp<Lielab::domain::sp>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexp", &Lielab::functions::dexp<Lielab::domain::su>, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &Lielab::functions::dexpinv<Lielab::domain::gl>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &Lielab::functions::dexpinv<Lielab::domain::rn>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &Lielab::functions::dexpinv<Lielab::domain::se>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &Lielab::functions::dexpinv<Lielab::domain::so>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &Lielab::functions::dexpinv<Lielab::domain::sp>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("dexpinv", &Lielab::functions::dexpinv<Lielab::domain::su>, "The dexpinv function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
    m_functions.def("left_product", &Lielab::functions::left_product, "Default action by left product.");

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
    m_transform.def("eanglebody123_to_dcm", &Lielab::transform::eanglebody123_to_dcm<double>, "The eanglebody123_to_dcm function.");
    m_transform.def("eanglebody231_to_dcm", &Lielab::transform::eanglebody231_to_dcm<double>, "The eanglebody231_to_dcm function.");
    m_transform.def("eanglebody312_to_dcm", &Lielab::transform::eanglebody312_to_dcm<double>, "The eanglebody312_to_dcm function.");
    m_transform.def("eanglebody132_to_dcm", &Lielab::transform::eanglebody132_to_dcm<double>, "The eanglebody132_to_dcm function.");
    m_transform.def("eanglebody213_to_dcm", &Lielab::transform::eanglebody213_to_dcm<double>, "The eanglebody213_to_dcm function.");
    m_transform.def("eanglebody321_to_dcm", &Lielab::transform::eanglebody321_to_dcm<double>, "The eanglebody321_to_dcm function.");
    m_transform.def("eanglebody121_to_dcm", &Lielab::transform::eanglebody121_to_dcm<double>, "The eanglebody121_to_dcm function.");
    m_transform.def("eanglebody131_to_dcm", &Lielab::transform::eanglebody131_to_dcm<double>, "The eanglebody131_to_dcm function.");
    m_transform.def("eanglebody212_to_dcm", &Lielab::transform::eanglebody212_to_dcm<double>, "The eanglebody212_to_dcm function.");
    m_transform.def("eanglebody232_to_dcm", &Lielab::transform::eanglebody232_to_dcm<double>, "The eanglebody232_to_dcm function.");
    m_transform.def("eanglebody313_to_dcm", &Lielab::transform::eanglebody313_to_dcm<double>, "The eanglebody313_to_dcm function.");
    m_transform.def("eanglebody323_to_dcm", &Lielab::transform::eanglebody323_to_dcm<double>, "The eanglebody323_to_dcm function.");
    m_transform.def("eanglespace123_to_dcm", &Lielab::transform::eanglespace123_to_dcm<double>, "The eanglespace123_to_dcm function.");
    m_transform.def("eanglespace231_to_dcm", &Lielab::transform::eanglespace231_to_dcm<double>, "The eanglespace231_to_dcm function.");
    m_transform.def("eanglespace312_to_dcm", &Lielab::transform::eanglespace312_to_dcm<double>, "The eanglespace312_to_dcm function.");
    m_transform.def("eanglespace132_to_dcm", &Lielab::transform::eanglespace132_to_dcm<double>, "The eanglespace132_to_dcm function.");
    m_transform.def("eanglespace213_to_dcm", &Lielab::transform::eanglespace213_to_dcm<double>, "The eanglespace213_to_dcm function.");
    m_transform.def("eanglespace321_to_dcm", &Lielab::transform::eanglespace321_to_dcm<double>, "The eanglespace321_to_dcm function.");
    m_transform.def("eanglespace121_to_dcm", &Lielab::transform::eanglespace121_to_dcm<double>, "The eanglespace121_to_dcm function.");
    m_transform.def("eanglespace131_to_dcm", &Lielab::transform::eanglespace131_to_dcm<double>, "The eanglespace131_to_dcm function.");
    m_transform.def("eanglespace212_to_dcm", &Lielab::transform::eanglespace212_to_dcm<double>, "The eanglespace212_to_dcm function.");
    m_transform.def("eanglespace232_to_dcm", &Lielab::transform::eanglespace232_to_dcm<double>, "The eanglespace232_to_dcm function.");
    m_transform.def("eanglespace313_to_dcm", &Lielab::transform::eanglespace313_to_dcm<double>, "The eanglespace313_to_dcm function.");
    m_transform.def("eanglespace323_to_dcm", &Lielab::transform::eanglespace323_to_dcm<double>, "The eanglespace323_to_dcm function.");
    m_transform.def("quaternion_to_dcm", &Lielab::transform::quaternion_to_dcm, "The quaternion_to_dcm function.");

    py::module m_topos = m.def_submodule("topos", "The topos submodule.");

    m_topos.def("Ad", &Lielab::topos::Ad);
    m_topos.def("dexp", &Lielab::topos::dexp, "The dexp function.", py::arg("a"), py::arg("b"), py::arg("order") = 5);
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


    // /*!
    // * Begin content for the "transform" submodule.
    // */
    // py::module m_transform = m.def_submodule("transform", "The transform submodule.");
    // m_transform.def("dcm_to_quaternion", &Lielab::transform::dcm_to_quaternion, "The dcm to quaternion transform function.");
    // m_transform.def("eanglebody123_to_dcm", &Lielab::transform::eanglebody123_to_dcm, "Euler angle body-123 rotation.");
    // m_transform.def("eanglebody231_to_dcm", &Lielab::transform::eanglebody231_to_dcm, "Euler angle body-231 rotation.");
    // m_transform.def("eanglebody312_to_dcm", &Lielab::transform::eanglebody312_to_dcm, "Euler angle body-312 rotation.");
    // m_transform.def("eanglebody132_to_dcm", &Lielab::transform::eanglebody132_to_dcm, "Euler angle body-132 rotation.");
    // m_transform.def("eanglebody213_to_dcm", &Lielab::transform::eanglebody213_to_dcm, "Euler angle body-213 rotation.");
    // m_transform.def("eanglebody321_to_dcm", &Lielab::transform::eanglebody321_to_dcm, "Euler angle body-321 rotation.");
    // m_transform.def("eanglebody121_to_dcm", &Lielab::transform::eanglebody121_to_dcm, "Euler angle body-121 rotation.");
    // m_transform.def("eanglebody131_to_dcm", &Lielab::transform::eanglebody131_to_dcm, "Euler angle body-131 rotation.");
    // m_transform.def("eanglebody212_to_dcm", &Lielab::transform::eanglebody212_to_dcm, "Euler angle body-212 rotation.");
    // m_transform.def("eanglebody232_to_dcm", &Lielab::transform::eanglebody232_to_dcm, "Euler angle body-232 rotation.");
    // m_transform.def("eanglebody313_to_dcm", &Lielab::transform::eanglebody313_to_dcm, "Euler angle body-313 rotation.");
    // m_transform.def("eanglebody323_to_dcm", &Lielab::transform::eanglebody323_to_dcm, "Euler angle body-323 rotation.");
    // m_transform.def("eanglespace123_to_dcm", &Lielab::transform::eanglespace123_to_dcm, "Euler angle space-123 rotation.");
    // m_transform.def("eanglespace231_to_dcm", &Lielab::transform::eanglespace231_to_dcm, "Euler angle space-231 rotation.");
    // m_transform.def("eanglespace312_to_dcm", &Lielab::transform::eanglespace312_to_dcm, "Euler angle space-312 rotation.");
    // m_transform.def("eanglespace132_to_dcm", &Lielab::transform::eanglespace132_to_dcm, "Euler angle space-132 rotation.");
    // m_transform.def("eanglespace213_to_dcm", &Lielab::transform::eanglespace213_to_dcm, "Euler angle space-213 rotation.");
    // m_transform.def("eanglespace321_to_dcm", &Lielab::transform::eanglespace321_to_dcm, "Euler angle space-321 rotation.");
    // m_transform.def("eanglespace121_to_dcm", &Lielab::transform::eanglespace121_to_dcm, "Euler angle space-121 rotation.");
    // m_transform.def("eanglespace131_to_dcm", &Lielab::transform::eanglespace131_to_dcm, "Euler angle space-131 rotation.");
    // m_transform.def("eanglespace212_to_dcm", &Lielab::transform::eanglespace212_to_dcm, "Euler angle space-212 rotation.");
    // m_transform.def("eanglespace232_to_dcm", &Lielab::transform::eanglespace232_to_dcm, "Euler angle space-232 rotation.");
    // m_transform.def("eanglespace313_to_dcm", &Lielab::transform::eanglespace313_to_dcm, "Euler angle space-313 rotation.");
    // m_transform.def("eanglespace323_to_dcm", &Lielab::transform::eanglespace323_to_dcm, "Euler angle space-323 rotation.");
    // m_transform.def("gibbs_to_dcm", &Lielab::transform::gibbs_to_dcm, "Gibbs to DCM.");
    // m_transform.def("dcm_to_gibbs", &Lielab::transform::dcm_to_gibbs, "DCM to Gibbs.");
    // m_transform.def("quaternion_to_dcm", &Lielab::transform::quaternion_to_dcm, "The quaternion to dcm transform function.");
}
