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

#include "pybind_CompositeGroup.hpp"

namespace py = pybind11;

void bind_CompositeGroup(py::module& m_domain)
{
    auto Lielab_domain_CompositeGroup = py::class_<Lielab::domain::CompositeGroup>(m_domain, "CompositeGroup");
    Lielab_domain_CompositeGroup.def_readwrite("point", &Lielab::domain::CompositeGroup::point);
    Lielab_domain_CompositeGroup.def(py::init());
    Lielab_domain_CompositeGroup.def(py::init<const Eigen::MatrixXcd&>());
    Lielab_domain_CompositeGroup.def_static("identity", &Lielab::domain::CompositeGroup::identity);
    Lielab_domain_CompositeGroup.def_static("project", &Lielab::domain::CompositeGroup::project);
    Lielab_domain_CompositeGroup.def(py::init<const int>());
    Lielab_domain_CompositeGroup.def(py::init<const std::vector<Lielab::domain::CompositeGroup::TYPES>&>());
    Lielab_domain_CompositeGroup.def("to_string", &Lielab::domain::CompositeGroup::to_string);
    Lielab_domain_CompositeGroup.def("get_dimension", &Lielab::domain::CompositeGroup::get_dimension);
    Lielab_domain_CompositeGroup.def("get_size", &Lielab::domain::CompositeGroup::get_size);
    Lielab_domain_CompositeGroup.def("is_abelian", &Lielab::domain::CompositeGroup::is_abelian);
    Lielab_domain_CompositeGroup.def("get_shape", &Lielab::domain::CompositeGroup::get_shape);
    Lielab_domain_CompositeGroup.def("__iter__",
        [](Lielab::domain::CompositeGroup& self)
        {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());
    Lielab_domain_CompositeGroup.def("get_dimensions", &Lielab::domain::CompositeGroup::get_dimensions);
    Lielab_domain_CompositeGroup.def("get_sizes", &Lielab::domain::CompositeGroup::get_sizes);
    Lielab_domain_CompositeGroup.def("get_shapes", &Lielab::domain::CompositeGroup::get_shapes);
    Lielab_domain_CompositeGroup.def("get_point", &Lielab::domain::CompositeGroup::get_point);
    Lielab_domain_CompositeGroup.def("serialize", &Lielab::domain::CompositeGroup::serialize);
    Lielab_domain_CompositeGroup.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::CompositeGroup::unserialize));
    Lielab_domain_CompositeGroup.def("get_matrix", &Lielab::domain::CompositeGroup::get_matrix);
    Lielab_domain_CompositeGroup.def("__call__",
        [](const Lielab::domain::CompositeGroup& self, const int index1, const int index2)
        {
            return self(index1, index2);
        });
    Lielab_domain_CompositeGroup.def("__getitem__",
        [](const Lielab::domain::CompositeGroup& self, const int index)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeGroup of length " + std::to_string(len));
            
            return self.point[_index];
        });
    Lielab_domain_CompositeGroup.def("__getitem__",
        [](const Lielab::domain::CompositeGroup& self, const py::slice slice)
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
                result.append(self.point[start + ii*step]);
            }
            return result;
        });
    Lielab_domain_CompositeGroup.def("__getitem__",
        [](const Lielab::domain::CompositeGroup& self, const py::slice slice)
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
                result.append(self.point[start + ii*step]);
            }
            return result;
        });
    Lielab_domain_CompositeGroup.def("__setitem__",
        [](Lielab::domain::CompositeGroup& self, const py::slice slice, const py::iterable& values)
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
                lielab_assert((index >= 0) && (index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeGroup of length " + std::to_string(len));
                indices[ii] = index;
            }

            auto it = values.begin();
            for (size_t ii = 0; ii < slicelength; ii++, it++)
            {
                self.point[indices[ii]] = (*it).cast<Lielab::domain::CompositeGroup::TYPES>();
            }
        });
    Lielab_domain_CompositeGroup.def("__delitem__",
        [](Lielab::domain::CompositeGroup& self, const int index)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeGroup of length " + std::to_string(len));
            
            self.point.erase(self.point.begin() + _index);
        });
    Lielab_domain_CompositeGroup.def("__delitem__",
        [](Lielab::domain::CompositeGroup& self, const py::slice slice)
        {
            const int len = static_cast<int>(self.point.size());

            size_t start, stop, step, slicelength;
            if (!slice.compute(len, &start, &stop, &step, &slicelength))
            {
                throw py::error_already_set();
            }

            std::vector<size_t> indices(slicelength);
            for (size_t ii = 0; ii < slicelength; ii++)
            {
                // Error check for out of bounds
                const size_t index = start + ii*step;
                lielab_assert((index >= 0) && (index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeGroup of length " + std::to_string(len));
                indices[ii] = index;
            }

            std::sort(indices.begin(), indices.end(), std::greater<size_t>());

            for (size_t index : indices)
            {
                self.point.erase(self.point.begin() + index);
            }
        });
    Lielab_domain_CompositeGroup.def(py::self * py::self);
    Lielab_domain_CompositeGroup.def(py::self *= py::self);
    Lielab_domain_CompositeGroup.def("inverse", &Lielab::domain::CompositeGroup::inverse);
    
    // other misc python
    Lielab_domain_CompositeGroup.def("__len__",
        [](const Lielab::domain::CompositeGroup& self)
        {
            return self.point.size();
        });
    Lielab_domain_CompositeGroup.def("append",
        [](Lielab::domain::CompositeGroup& self, const py::object& value)
        {
            self.point.push_back(value.cast<Lielab::domain::CompositeGroup::TYPES>());
        });
    Lielab_domain_CompositeGroup.def("extend",
        [](Lielab::domain::CompositeGroup& self, const py::iterable& values)
        {
            const size_t valueslength = std::distance(values.begin(), values.end());
            auto it = values.begin();
            for (size_t ii = 0; ii < valueslength; ii++, it++)
            {
                self.point.push_back((*it).cast<Lielab::domain::CompositeGroup::TYPES>());
            }
        });
    // TODO: Leave this off for now. CompositeAlgebra logically collides with addition and this might reinforce bad habits.
    // Lielab_domain_CompositeGroup.def("__iadd__",
    //     [](Lielab::domain::CompositeGroup& self, const py::iterable& values)
    //     {
    //         const size_t valueslength = std::distance(values.begin(), values.end());
    //         auto it = values.begin();
    //         for (size_t ii = 0; ii < valueslength; ii++, it++)
    //         {
    //             self.point.push_back((*it).cast<Lielab::domain::CompositeGroup::TYPES>());
    //         }
    //         return self;
    //     });
    Lielab_domain_CompositeGroup.def("insert",
        [](Lielab::domain::CompositeGroup& self, const int index, const py::object& value)
        {
            self.point.insert(self.point.begin() + index, value.cast<Lielab::domain::CompositeGroup::TYPES>());
        });
    Lielab_domain_CompositeGroup.def("pop",
        [](Lielab::domain::CompositeGroup& self, const int index)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeGroup of length " + std::to_string(len));

            const auto out = self.point[_index];
            self.point.erase(self.point.begin() + _index);
            return out;
        });
    Lielab_domain_CompositeGroup.def("clear",
        [](Lielab::domain::CompositeGroup& self)
        {
            self.point.clear();
        });
    Lielab_domain_CompositeGroup.def(py::pickle([](const Lielab::domain::CompositeGroup& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point);
        },
        [](py::tuple t)
        {
            // __setstate__
            lielab_assert(t.size() == 1, "Invalid state for CompositeGroup.");
            Lielab::domain::CompositeGroup obj;
            obj.point = t[0].cast<std::vector<Lielab::domain::CompositeGroup::TYPES>>();
            return obj;
        }));
    Lielab_domain_CompositeGroup.def("__repr__",
        [](const Lielab::domain::CompositeGroup & self)
        {
            return "<lielab.domain.CompositeGroup>";
        });
    Lielab_domain_CompositeGroup.def("__str__",
        [](const Lielab::domain::CompositeGroup & self)
        {
            return "<lielab.domain.CompositeGroup>";
        });
}
