#include "CompositeManifold.hpp"

#include "../liealgebras.hpp"
#include "../liegroups.hpp"

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

CompositeManifold::CompositeManifold()
{
    
}

// CompositeManifold::~CompositeManifold();

CompositeManifold::CompositeManifold(const int n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{CompositeManifold} \f}
    *
    * Constructor instantiating a \f$\mathfrak{CompositeManifold}\f$ as multiple empty glc objects.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::CompositeManifold x(3), y(4), z(5);
    * 
    * @param[in] n
    */

    this->point.reserve(n);
    for (int ii = 0; ii < n; ii++)
    {
        this->point.push_back(glc(0));
    }
}

CompositeManifold::CompositeManifold(std::initializer_list<CompositeManifold::TYPES> others)
{
    /*!
    * Initializer list constructor for CompositeManifold
    *
    * Enables construction like:
    *     Lielab::domain::RN R(3);
    *     Lielab::domain::SO O(3);
    *     Lielab::domain::CompositeManifold M{R, O};
    */

    this->point = std::vector<CompositeManifold::TYPES>{std::move(others)};
}

CompositeManifold::CompositeManifold(const std::vector<CompositeManifold::TYPES>& others)
{
    /*!
    * Vector list constructor for CompositeManifold
    *
    * Not needed for C++, but enables construction in Python like:
    *     R = lielab.domain.RN(3)
    *     O = lielab.domain.SO(3)
    *     M = lielab.domain.CompositeManifold([R, O])
    */

    this->point = others;
}

std::string CompositeManifold::to_string() const
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

        if (ii < (sz - 1)) out += " x ";
        ii++;
    }

    return out;
}

int CompositeManifold::get_dimension() const
{
    const std::vector<int> dims = this->get_dimensions();
    return std::accumulate(dims.begin(), dims.end(), 0);
}

int CompositeManifold::get_size() const
{
    const std::vector<int> sizes = this->get_sizes();
    return std::accumulate(sizes.begin(), sizes.end(), 0);
}

std::vector<CompositeManifold::TYPES>::iterator CompositeManifold::begin()
{
    return this->point.begin();
}

std::vector<CompositeManifold::TYPES>::iterator CompositeManifold::end()
{
    return this->point.end();
}

std::vector<CompositeManifold::TYPES>::const_iterator CompositeManifold::begin() const
{
    return this->point.begin();
}

std::vector<CompositeManifold::TYPES>::const_iterator CompositeManifold::end() const
{
    return this->point.end();
}

std::vector<int> CompositeManifold::get_dimensions() const
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

std::vector<int> CompositeManifold::get_sizes() const
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

CompositeManifold::point_t CompositeManifold::get_point() const
{
    return this->point;
}

Eigen::VectorXd CompositeManifold::serialize() const
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

void CompositeManifold::unserialize(const Eigen::VectorXd& vector)
{
    int start = 0;
    for (auto& element : this->point)
    {
        std::visit([&](auto& _element)
        {
            const int size = _element.get_size();
            _element.unserialize(vector(Eigen::seqN(start, size)));
            start += size;
        }, element);
    }
}

void CompositeManifold::unserialize(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->unserialize(Eigen::VectorXd{std::move(vector)});
}

const CompositeManifold::dataproxy CompositeManifold::operator[](const int index) const
{
    const int len = static_cast<int>(this->point.size());

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? len + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeManifold of length " + std::to_string(len));

    return CompositeManifold::dataproxy{const_cast<CompositeManifold::TYPES&>(this->point[_index])};
}

CompositeManifold::dataproxy CompositeManifold::operator[](const int index)
{
    const int len = static_cast<int>(this->point.size());

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? len + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < len), "Index " + std::to_string(index) + " is out of bounds for CompositeManifold of length " + std::to_string(len));

    return CompositeManifold::dataproxy{this->point[_index]};
}

}
