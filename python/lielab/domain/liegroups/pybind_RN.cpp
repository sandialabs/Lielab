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

#include "pybind_RN.hpp"

namespace py = pybind11;

void bind_RN(py::module& m_domain)
{
    auto Lielab_domain_RN = py::class_<Lielab::domain::RN>(m_domain, "RN");
    Lielab_domain_RN.def_readwrite("point", &Lielab::domain::RN::point);
    Lielab_domain_RN.def(py::init<>());
    Lielab_domain_RN.def(py::init<const Eigen::MatrixXd&>());
    Lielab_domain_RN.def_static("identity", &Lielab::domain::RN::identity);
    Lielab_domain_RN.def_static("project", &Lielab::domain::RN::project);
    Lielab_domain_RN.def(py::init<const int>());
    Lielab_domain_RN.def_static("from_vector", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::RN::from_vector));
    Lielab_domain_RN.def("to_string", &Lielab::domain::RN::to_string);
    Lielab_domain_RN.def("__repr__", [](const Lielab::domain::RN& self)
        {
            return "<lielab.domain.RN>";
        });
    Lielab_domain_RN.def("__str__", [](const Lielab::domain::RN& self)
        {
            return "<lielab.domain.RN>";
        });
    Lielab_domain_RN.def("get_dimension", &Lielab::domain::RN::get_dimension);
    Lielab_domain_RN.def("get_size", &Lielab::domain::RN::get_size);
    Lielab_domain_RN.def("is_abelian", &Lielab::domain::RN::is_abelian);
    Lielab_domain_RN.def("get_shape", &Lielab::domain::RN::get_shape);
    Lielab_domain_RN.def("get_point", &Lielab::domain::RN::get_point);
    Lielab_domain_RN.def("serialize", &Lielab::domain::RN::serialize);
    Lielab_domain_RN.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::RN::unserialize));
    Lielab_domain_RN.def("get_matrix", &Lielab::domain::RN::get_matrix);
    Lielab_domain_RN.def("__call__", [](const Lielab::domain::RN& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_RN.def("__getitem__",
        [](const Lielab::domain::RN& self, const int index)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for RN of length " + std::to_string(len));
            
            return self[index];
        });
    Lielab_domain_RN.def("__getitem__",
        [](const Lielab::domain::RN& self, const py::slice slice)
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
    Lielab_domain_RN.def("__setitem__",
        [](Lielab::domain::RN& self, const int index, const py::object& value)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for RN of length " + std::to_string(len));

            self[_index] = value.cast<Lielab::domain::RN::field_t>();
        });
    Lielab_domain_RN.def("__setitem__",
        [](Lielab::domain::RN& self, const py::slice slice, const py::iterable& values)
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
                lielab_assert((index >= 0) && (index < len), "Index " + std::to_string(index) + " is out of bounds for RN of length " + std::to_string(len));
                indices[ii] = index;
            }

            auto it = values.begin();
            for (size_t ii = 0; ii < slicelength; ii++, it++)
            {
                self[indices[ii]] = (*it).cast<Lielab::domain::RN::field_t>();
            }
        });
    Lielab_domain_RN.def(py::self * py::self);
    Lielab_domain_RN.def(py::self *= py::self);
    Lielab_domain_RN.def("inverse", &Lielab::domain::RN::inverse);

    // Misc python functions
    Lielab_domain_RN.def(py::pickle([](const Lielab::domain::RN& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point, obj._shape);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 2)
                throw std::runtime_error("RN: Invalid state.");

            Lielab::domain::RN obj;
            obj.point = t[0].cast<Lielab::domain::RN::point_t>();
            obj._shape = t[1].cast<int>();
            return obj;
        }));
}
