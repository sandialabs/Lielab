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

#include "pybind_Grassmannian.hpp"

namespace py = pybind11;

void bind_Grassmannian(py::module& m_domain)
{
    auto Lielab_domain_Grassmannian = py::class_<Lielab::domain::Grassmannian>(m_domain, "Grassmannian");
    Lielab_domain_Grassmannian.def_readwrite("point", &Lielab::domain::Grassmannian::point);
    Lielab_domain_Grassmannian.def_readwrite("axes", &Lielab::domain::Grassmannian::axes);
    Lielab_domain_Grassmannian.def_readwrite("point_projector", &Lielab::domain::Grassmannian::point_projector);
    Lielab_domain_Grassmannian.def(py::init());
    Lielab_domain_Grassmannian.def(py::init<const int, const int>());
    Lielab_domain_Grassmannian.def(py::init<const Eigen::VectorXd&, const Eigen::MatrixXd&>());
    Lielab_domain_Grassmannian.def_static("project", &Lielab::domain::Grassmannian::project);
    Lielab_domain_Grassmannian.def("to_string", &Lielab::domain::Grassmannian::to_string);
    Lielab_domain_Grassmannian.def("__repr__",
        [](const Lielab::domain::Grassmannian& self)
        {
            return "<lielab.domain.Grassmannian>";
        });
    Lielab_domain_Grassmannian.def("__str__",
        [](const Lielab::domain::Grassmannian& self)
        {
            return "<lielab.domain.Grassmannian>";
        });
    Lielab_domain_Grassmannian.def("get_dimension", &Lielab::domain::Grassmannian::get_dimension);
    Lielab_domain_Grassmannian.def("get_size", &Lielab::domain::Grassmannian::get_size);
    Lielab_domain_Grassmannian.def("get_point", &Lielab::domain::Grassmannian::get_point);
    Lielab_domain_Grassmannian.def("serialize", &Lielab::domain::Grassmannian::serialize);
    Lielab_domain_Grassmannian.def("unserialize", py::overload_cast<const Eigen::VectorXd&>(&Lielab::domain::Grassmannian::unserialize));
    Lielab_domain_Grassmannian.def("__getitem__",
        [](const Lielab::domain::Grassmannian& self, const int index)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for Grassmannian of size " + std::to_string(len));
            
            return self.point(_index);
        });
    Lielab_domain_Grassmannian.def("__getitem__",
        [](const Lielab::domain::Grassmannian& self, const py::slice slice)
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
                result.append(self.point(start + ii*step));
            }
            return result;
        });
    Lielab_domain_Grassmannian.def("__setitem__",
        [](Lielab::domain::Grassmannian& self, const int index, const py::object& value)
        {
            const int len = static_cast<int>(self.point.size());

            // If input index is negative, index from the back of the array
            const int _index = (index < 0) ? len + index : index;

            // Error check for out of bounds
            lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for Grassmannian of size " + std::to_string(len));

            self.point(_index) = value.cast<double>();
        });
    Lielab_domain_Grassmannian.def("__setitem__",
        [](Lielab::domain::Grassmannian& self, const py::slice slice, const py::iterable& values)
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
                lielab_assert((index >= 0) && (index < len), "Index " + std::to_string(index) + " is out of bounds for Grassmannian of size " + std::to_string(len));
                indices[ii] = index;
            }

            auto it = values.begin();
            for (size_t ii = 0; ii < slicelength; ii++, it++)
            {
                self.point(indices[ii]) = (*it).cast<double>();
            }
        });
    Lielab_domain_Grassmannian.def("project_point", &Lielab::domain::Grassmannian::project_point);
    Lielab_domain_Grassmannian.def("project_vector_onto_tangent_space", &Lielab::domain::Grassmannian::project_vector_onto_tangent_space);
    Lielab_domain_Grassmannian.def("project_vector_onto_normal_space", &Lielab::domain::Grassmannian::project_vector_onto_normal_space);
    Lielab_domain_Grassmannian.def("axes_intersection", &Lielab::domain::Grassmannian::axes_intersection);
    
    // Other misc python methods
    Lielab_domain_Grassmannian.def(py::pickle([](const Lielab::domain::Grassmannian& obj)
        {
            // __getstate__
            return py::make_tuple(obj.point, obj.axes);
        },
        [](py::tuple t)
        {
            // __setstate__
            if (t.size() != 2)
                throw std::runtime_error("Grassmannian: Invalid state.");

            Lielab::domain::Grassmannian obj;
            obj.point = t[0].cast<Lielab::domain::Grassmannian::point_t>();
            obj.axes = t[1].cast<Eigen::MatrixXd>();
            return obj;
        }));
}
