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

#include "pybind_CN.hpp"

namespace py = pybind11;

void bind_CN(py::module& m_domain)
{
    auto Lielab_domain_CN = py::class_<Lielab::domain::CN>(m_domain, "CN");
    Lielab_domain_CN.def_readwrite("point", &Lielab::domain::CN::point);
    Lielab_domain_CN.def(py::init<>());
    Lielab_domain_CN.def(py::init<const Eigen::MatrixXcd&>());
    Lielab_domain_CN.def_static("identity", &Lielab::domain::CN::identity);
    Lielab_domain_CN.def_static("project", &Lielab::domain::CN::project);
    Lielab_domain_CN.def(py::init<const int>());
    Lielab_domain_CN.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::CN::from_vector));
    Lielab_domain_CN.def_static("from_complex_vector", py::overload_cast<const Eigen::VectorXcd&>(&Lielab::domain::CN::from_complex_vector));
    Lielab_domain_CN.def("to_string", &Lielab::domain::CN::to_string);
    Lielab_domain_CN.def("get_dimension", &Lielab::domain::CN::get_dimension);
    Lielab_domain_CN.def("get_size", &Lielab::domain::CN::get_size);
    Lielab_domain_CN.def("is_abelian", &Lielab::domain::CN::is_abelian);
    Lielab_domain_CN.def("get_shape", &Lielab::domain::CN::get_shape);
    Lielab_domain_CN.def("get_point", &Lielab::domain::CN::get_point);
    Lielab_domain_CN.def("serialize", &Lielab::domain::CN::serialize);
    Lielab_domain_CN.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::CN::unserialize));
    Lielab_domain_CN.def("get_matrix", &Lielab::domain::CN::get_matrix);
    Lielab_domain_CN.def("__call__", [](const Lielab::domain::CN& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_CN.def("to_complex_vector", &Lielab::domain::CN::to_complex_vector);
    Lielab_domain_CN.def("__call__", [](const Lielab::domain::CN& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_CN.def("__getitem__",
        [](const Lielab::domain::CN& self, const int index)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CN of length " + std::to_string(len));
            
            return self[index];
        });
    Lielab_domain_CN.def("__getitem__",
        [](const Lielab::domain::CN& self, const py::slice slice)
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
    Lielab_domain_CN.def("__setitem__",
        [](Lielab::domain::CN& self, const int index, const py::object& value)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CN of length " + std::to_string(len));

            self[_index] = value.cast<Lielab::domain::CN::field_t>();
        });
    Lielab_domain_CN.def("__setitem__",
        [](Lielab::domain::CN& self, const py::slice slice, const py::iterable& values)
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
                lielab_assert((index >= 0) && (index < len), "Index " + std::to_string(index) + " is out of bounds for CN of length " + std::to_string(len));
                indices[ii] = index;
            }

            auto it = values.begin();
            for (size_t ii = 0; ii < slicelength; ii++, it++)
            {
                self[indices[ii]] = (*it).cast<Lielab::domain::CN::field_t>();
            }
        });
    Lielab_domain_CN.def(py::self * py::self);
    Lielab_domain_CN.def(py::self *= py::self);
    Lielab_domain_CN.def("inverse", &Lielab::domain::CN::inverse);

    // Extra python funcs
    Lielab_domain_CN.def(py::pickle([](const Lielab::domain::CN& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point, obj._shape);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 2)
                throw std::runtime_error("CN: Invalid state.");

            Lielab::domain::CN obj;
            obj.point = t[0].cast<Lielab::domain::CN::point_t>();
            obj._shape = t[1].cast<int>();
            return obj;
        }));
    Lielab_domain_CN.def("__repr__", [](const Lielab::domain::CN& self)
        {
            return "<lielab.domain.CN>";
        });
    Lielab_domain_CN.def("__str__", [](const Lielab::domain::CN& self)
        {
            return "<lielab.domain.CN>";
        });
}
