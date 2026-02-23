#include "CompositeGroup.hpp"

#include "Lielab/domain/liealgebras.hpp" // TODO: Why do we need the algebras? Remove this
#include "Lielab/domain/liegroups.hpp"

#include "Lielab/testing.hpp"
#include "Lielab/utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <variant>

namespace Lielab::domain
{

CompositeGroup::CompositeGroup()
{
    
}

// CompositeGroup::~CompositeGroup();

CompositeGroup::CompositeGroup(const CompositeGroup::matrix_t& matrix)
{
    this->point.push_back(GLC(matrix));
}

CompositeGroup CompositeGroup::identity(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{CompositeGroup} \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The CompositeGroup element. 
    */

    return CompositeGroup({GLC::identity(shape)});
}

CompositeGroup CompositeGroup::project(const CompositeGroup::matrix_t& matrix)
{
    return CompositeGroup({GLC::project(matrix)});
}

CompositeGroup::CompositeGroup(const int n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{CompositeGroup} \f}
    *
    * Constructor instantiating a \f$\mathfrak{CompositeGroup}\f$ as multiple empty GLC objects.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::CompositeGroup x(3), y(4), z(5);
    * 
    * @param[in] n
    */

    this->point.reserve(n);
    for (int ii = 0; ii < n; ii++)
    {
        this->point.push_back(GLC(0));
    }
}

CompositeGroup::CompositeGroup(std::initializer_list<TYPES> others)
{
    /*!
    * Initializer list constructor for CompositeGroup
    *
    * Enables construction like:
    *     Lielab::domain::RN R(3);
    *     Lielab::domain::SO O(3);
    *     Lielab::domain::CompositeGroup M{R, O};
    */

    this->point = std::vector<TYPES>{std::move(others)};
}

CompositeGroup::CompositeGroup(const std::vector<TYPES>& others)
{
    /*!
    * Vector list constructor for CompositeGroup
    *
    * Not needed for C++, but enables construction in Python like:
    *     R = lielab.domain.RN(3)
    *     O = lielab.domain.SO(3)
    *     M = lielab.domain.CompositeGroup([R, O])
    */

    this->point = others;
}

std::string CompositeGroup::to_string() const
{
    std::string out = "";
    const int sz = static_cast<int>(this->point.size());

    int ii = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            out += _element.to_string();
        }, element);

        if (ii < sz - 1) out += " x ";
        ii++;
    }

    return out;
}

int CompositeGroup::get_dimension() const
{
    const std::vector<int> dims = this->get_dimensions();
    return std::accumulate(dims.begin(), dims.end(), 0);
}

int CompositeGroup::get_size() const
{
    const std::vector<int> sizes = this->get_sizes();
    return std::accumulate(sizes.begin(), sizes.end(), 0);
}

bool CompositeGroup::is_abelian() const
{
    bool abelian = true;
    for (auto& element : this->point)
    {
        std::visit([&](auto& _element)
        {
            if (!_element.is_abelian()) abelian = false;
        }, element);
    }
    return abelian;
}

int CompositeGroup::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    const std::vector<int> shapes = this->get_shapes();
    return std::accumulate(shapes.begin(), shapes.end(), 0);
}

std::vector<CompositeGroup::TYPES>::iterator CompositeGroup::begin()
{
    return this->point.begin();
}

std::vector<CompositeGroup::TYPES>::iterator CompositeGroup::end()
{
    return this->point.end();
}

std::vector<CompositeGroup::TYPES>::const_iterator CompositeGroup::begin() const
{
    return this->point.begin();
}

std::vector<CompositeGroup::TYPES>::const_iterator CompositeGroup::end() const
{
    return this->point.end();
}

std::vector<int> CompositeGroup::get_dimensions() const
{
    std::vector<int> dimensions(this->point.size());

    int ii = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            dimensions[ii] = _element.get_dimension();
        }, element);
        ii++;
    }

    return dimensions;
}

std::vector<int> CompositeGroup::get_sizes() const
{
    std::vector<int> sizes(this->point.size());

    int ii = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            sizes[ii] = _element.get_size();
        }, element);
        ii++;
    }

    return sizes;
}

std::vector<int> CompositeGroup::get_shapes() const
{
    std::vector<int> shapes(this->point.size());

    int ii = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            shapes[ii] = _element.get_shape();
        }, element);
        ii++;
    }

    return shapes;
}

CompositeGroup::point_t CompositeGroup::get_point() const
{
    return this->point;
}

Eigen::VectorXd CompositeGroup::serialize() const
{
    const int length = static_cast<int>(this->point.size());

    if (length == 0) return Eigen::VectorXd(0);
    
    std::vector<Eigen::VectorXd> serials(length);
    
    int index = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            serials[index] = _element.serialize();
        }, element);
        index++;
    }
    
    Eigen::VectorXd out = Lielab::utils::concatenate(serials);// TODO: Remove this
    return out;
}

void CompositeGroup::unserialize(const Eigen::VectorXd& serialized)
{
    int start = 0;
    for (auto& element : this->point)
    {
        std::visit([&](auto& _element)
        {
            const int size = _element.get_size();
            _element.unserialize(serialized(Eigen::seqN(start, size)));
            start += size;
        }, element);
    }
}

void CompositeGroup::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

CompositeGroup::matrix_t CompositeGroup::get_matrix() const
{
    const int shape = this->get_shape();
    Eigen::MatrixXcd out = Eigen::MatrixXcd::Zero(shape, shape);
    
    int start = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            const int _shape = _element.get_shape();
            out(Eigen::seqN(start, _shape), Eigen::seqN(start, _shape)) = _element.get_matrix();
            start += _shape;
        }, element);
    }

    return out;
}

CompositeGroup::field_t CompositeGroup::operator()(const int index1, const int index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int shape = this->get_shape();
    if (shape == 0) return std::complex<double>(nan, nan);

    // If input index is negative, index from the back of the array
    const int _index1 = (index1 < 0) ? shape + index1 : index1;
    const int _index2 = (index2 < 0) ? shape + index2 : index2;

    // Error check for out of bounds
    if (_index1 < 0) return std::complex<double>(nan, nan);
    if (_index1 >= shape) return std::complex<double>(nan, nan);
    if (_index2 < 0) return std::complex<double>(nan, nan);
    if (_index2 >= shape) return std::complex<double>(nan, nan);
    
    const std::vector<int> shapes = this->get_shapes();

    int relind = 0;
    for (int ii = 0; ii < static_cast<int>(shapes.size()); ii++)
    {
        // Sparse components
        if ((_index1 - relind) < 0 || (_index2 - relind) < 0)
        {
            return std::complex<double>(0.0, 0.0);
        }

        if ((_index1 - relind) < shapes[ii] && (_index2 - relind) < shapes[ii])
        {
            return std::visit([&](const auto& _element)
            {
                return static_cast<std::complex<double>>(_element(_index1 - relind, _index2 - relind));
            }, this->point[ii]);
        }

        relind += shapes[ii];
    }

    return std::complex<double>(nan, nan);
}

const CompositeGroup::dataproxy CompositeGroup::operator[](const int index) const
{
    const int len = static_cast<int>(this->point.size());

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? len + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeGroup of length " + std::to_string(len));

    return CompositeGroup::dataproxy{const_cast<CompositeGroup::TYPES&>(this->point[_index])};
}

CompositeGroup::dataproxy CompositeGroup::operator[](const int index)
{
    const int len = static_cast<int>(this->point.size());

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? len + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeGroup of length " + std::to_string(len));

    return CompositeGroup::dataproxy{this->point[_index]};
}

CompositeGroup CompositeGroup::operator*(const CompositeGroup& other) const
{
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(*this, other), "Unable to take product of topologically inconsistent groups: (" + this->to_string() + ") !≅ (" + other.to_string() + ").");

    CompositeGroup out;
    int index = 0;
    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            out.point.push_back(_element*std::get<other_t>(other.point[index]));
        }, element);
        index++;
    }

    return out;
}

CompositeGroup& CompositeGroup::operator*=(const CompositeGroup& other)
{
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(*this, other), "Unable to take product of topologically inconsistent groups: (" + this->to_string() + ") !≅ (" + other.to_string() + ").");

    int index = 0;
    for (auto& element : this->point)
    {
        std::visit([&](auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            _element *= std::get<other_t>(other.point[index]);
        }, element);
        index++;
    }

    return *this;
}

CompositeGroup CompositeGroup::inverse() const
{
    CompositeGroup out;

    for (const auto& element : this->point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(_element.inverse());
        }, element);
    }

    return out;
}

}
