#ifndef LIELAB_DOMAIN_LIEGROUP_TPP
#define LIELAB_DOMAIN_LIEGROUP_TPP

#include "LieGroup.hpp"

#include <Eigen/Core>

#include <string>

namespace Lielab::domain
{

template <typename Field>
bool LieGroup<Field>::is_abelian() const
{
    return false;
}

template <typename Field>
bool LieGroup<Field>::is_complex() const
{
    return true;
}

template <typename Field>
std::string LieGroup<Field>::to_string() const
{
    if (this->is_complex())
    {
        return "G(" + std::to_string(this->get_shape()) + ", C)";
    }
    
    return "G(" + std::to_string(this->get_shape()) + ", R)";
}

template <typename Field>
LieGroup<Field>::LieGroup()
{

}

template <typename Field>
LieGroup<Field>::LieGroup(const size_t n)
{
    
}

template <typename Field>
LieGroup<Field>::LieGroup(const LieGroup::matrix_t& other)
{
    assert(other.rows() == other.cols());
    this->data = other;
}

template <typename Field>
LieGroup<Field>::~LieGroup()
{

}

}

#endif
