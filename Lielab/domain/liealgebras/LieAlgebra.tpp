#ifndef LIELAB_DOMAIN_LIEALGEBRAS_LIEALGEBRA_TPP
#define LIELAB_DOMAIN_LIEALGEBRAS_LIEALGEBRA_TPP

#include "LieAlgebra.hpp"

#include <Eigen/Core>

#include <string>

namespace Lielab::domain
{

template <typename Field>
bool LieAlgebra<Field>::is_abelian() const
{
    return false;
}

template <typename Field>
bool LieAlgebra<Field>::is_complex() const
{
    return true;
}

template <typename Field>
std::string LieAlgebra<Field>::to_string() const
{
    if (this->is_complex())
    {
        return "g(" + std::to_string(this->get_shape()) + ", C)";
    }
    
    return "g(" + std::to_string(this->get_shape()) + ", R)";
}

template <typename Field>
LieAlgebra<Field>::LieAlgebra()
{

}

template <typename Field>
LieAlgebra<Field>::LieAlgebra(const size_t n)
{
    
}

template <typename Field>
LieAlgebra<Field>::LieAlgebra(const LieAlgebra::matrix_t& other)
{
    assert(other.rows() == other.cols());
    this->data = other;
}

template <typename Field>
LieAlgebra<Field>::~LieAlgebra()
{

}

template <typename Field>
LieAlgebra<Field>::data_t LieAlgebra<Field>::get_data() const
{
    return this->data;
}

// template <typename Field>
// LieAlgebra<Field>::matrix_t LieAlgebra<Field>::get_matrix() const
// {
//     return this->data;
// }

// get_vector
// set_vector

template <typename Field>
size_t LieAlgebra<Field>::get_shape() const
{
    return this->data.rows();
}

}

#endif
