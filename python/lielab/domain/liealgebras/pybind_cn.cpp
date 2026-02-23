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

void bind_cn(py::module& m_domain)
{
    auto Lielab_domain_cn = py::class_<Lielab::domain::cn>(m_domain, "cn");
    Lielab_domain_cn.def_readwrite("point", &Lielab::domain::cn::point);
    Lielab_domain_cn.def(py::init<>());
    Lielab_domain_cn.def(py::init<const Eigen::MatrixXcd&>());
    Lielab_domain_cn.def_static("basis", &Lielab::domain::cn::basis);
    Lielab_domain_cn.def_static("zero", &Lielab::domain::cn::zero);
    Lielab_domain_cn.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::cn::from_vector));
    Lielab_domain_cn.def_static("project", &Lielab::domain::cn::project);
    Lielab_domain_cn.def(py::init<const int>());
    Lielab_domain_cn.def_static("from_complex_vector", py::overload_cast<const Eigen::VectorXcd&>(&Lielab::domain::cn::from_complex_vector));
    Lielab_domain_cn.def("to_string", &Lielab::domain::cn::to_string);
    Lielab_domain_cn.def("__repr__", [](const Lielab::domain::cn& self)
        {
            return "<lielab.domain.cn>";
        });
    Lielab_domain_cn.def("__str__", [](const Lielab::domain::cn& self)
        {
            return "<lielab.domain.cn>";
        });
    Lielab_domain_cn.def("get_dimension", &Lielab::domain::cn::get_dimension);
    Lielab_domain_cn.def("get_size", &Lielab::domain::cn::get_size);
    Lielab_domain_cn.def("is_abelian", &Lielab::domain::cn::is_abelian);
    Lielab_domain_cn.def("get_shape", &Lielab::domain::cn::get_shape);
    Lielab_domain_cn.def("get_point", &Lielab::domain::cn::get_point);
    Lielab_domain_cn.def("serialize", &Lielab::domain::cn::serialize);
    Lielab_domain_cn.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::cn::unserialize));
    Lielab_domain_cn.def("get_matrix", &Lielab::domain::cn::get_matrix);
    Lielab_domain_cn.def("get_vector", &Lielab::domain::cn::get_vector);
    Lielab_domain_cn.def("set_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::cn::set_vector));
    Lielab_domain_cn.def("__call__", [](const Lielab::domain::cn& self, const int index)
        {
            return self(index);
        });
    Lielab_domain_cn.def("__call__", [](const Lielab::domain::cn& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_cn.def("to_complex_vector", &Lielab::domain::cn::to_complex_vector);
    Lielab_domain_cn.def("__getitem__",
        [](const Lielab::domain::cn& self, const int index)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for cn of length " + std::to_string(len));
            
            return self[index];
        });
    Lielab_domain_cn.def("__getitem__",
        [](const Lielab::domain::cn& self, const py::slice slice)
        {
            const int len = static_cast<int>(self.point.size());

            size_t start, stop, step, slicelength;
            if (!slice.compute(len, &start, &stop, &step, &slicelength))
            {
                throw py::error_already_set();
            }

            py::list result;
            for (size_t ii = 0; ii < slicelength; ii++)
            {
                result.append(self[start + ii*step]);
            }
            return result;
        });
    Lielab_domain_cn.def("__setitem__",
        [](Lielab::domain::cn& self, const int index, const py::object& value)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for cn of length " + std::to_string(len));

            self[_index] = value.cast<Lielab::domain::cn::field_t>();
        });
    Lielab_domain_cn.def("__setitem__",
        [](Lielab::domain::cn& self, const py::slice slice, const py::iterable& values)
        {
            const int len = static_cast<int>(self.point.size());

            size_t start, stop, step, slicelength;
            if (!slice.compute(len, &start, &stop, &step, &slicelength))
            {
                throw py::error_already_set();
            }

            const size_t valueslength = std::distance(values.begin(), values.end());
            lielab_assert(slicelength == valueslength, "Slice length (" + std::to_string(slicelength) + ") does not match number of given values (" + std::to_string(valueslength) + ").");

            std::vector<size_t> indices(slicelength);
            for (size_t ii = 0; ii < slicelength; ii++)
            {
                // Error check for out of bounds
                const size_t index = start + ii*step;
                lielab_assert((index >= 0) && (index < len), "Index " + std::to_string(index) + " is out of bounds for cn of length " + std::to_string(len));
                indices[ii] = index;
            }

            auto it = values.begin();
            for (size_t ii = 0; ii < slicelength; ii++, it++)
            {
                self[indices[ii]] = (*it).cast<Lielab::domain::cn::field_t>();
            }
        });
    Lielab_domain_cn.def(py::self + py::self);
    Lielab_domain_cn.def(py::self += py::self);
    Lielab_domain_cn.def(py::self - py::self);
    Lielab_domain_cn.def(py::self -= py::self);
    Lielab_domain_cn.def(-py::self);
    Lielab_domain_cn.def(py::self * int());
    Lielab_domain_cn.def(py::self * double());
    Lielab_domain_cn.def(py::self *= int());
    Lielab_domain_cn.def(int() * py::self);
    Lielab_domain_cn.def(double() * py::self);
    Lielab_domain_cn.def(py::self *= double());
    Lielab_domain_cn.def(py::self / int());
    Lielab_domain_cn.def(py::self / double());
    Lielab_domain_cn.def(py::self /= int());
    Lielab_domain_cn.def(py::self /= double());
    Lielab_domain_cn.def(py::self * std::complex<double>());
    Lielab_domain_cn.def(std::complex<double>() * py::self);
    Lielab_domain_cn.def(py::self *= std::complex<double>());
    Lielab_domain_cn.def(py::self / std::complex<double>());
    Lielab_domain_cn.def(py::self /= std::complex<double>());

    // Other misc python
    Lielab_domain_cn.def(py::pickle([](const Lielab::domain::cn& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point, obj._shape);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 2)
                throw std::runtime_error("cn: Invalid state.");

            Lielab::domain::cn obj;
            obj.point = t[0].cast<Lielab::domain::cn::point_t>();
            obj._shape = t[1].cast<int>();
            return obj;
        }));
}
