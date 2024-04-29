#ifndef _LIELAB_FUNCTIONS_ACTIONS_HPP
#define _LIELAB_FUNCTIONS_ACTIONS_HPP

#include "exp.hpp"

namespace Lielab
{
namespace functions
{
Lielab::domain::CompositeManifold left_product(const Lielab::domain::CompositeManifold & Left, const Lielab::domain::CompositeManifold & Right)
{
    /*!
     * Default action by left product.
     *
     * TODO: Probably belongs in the topos namespace.
     *
     * @param[in] Left
     * @param[in] Right
     * @param[out] out Left*Right
     * 
     * Author / Date: Sparapany / 2022
     */

    // Simple error checking on the inputs
    // TODO: Move this error checking to CompositeManifold?
    if (Left.space.size() != Right.space.size())
    {
        throw SizeError("left_product: Spaces of CompositeManifold must be the same size (" + std::to_string(Left.space.size()) + " != " + std::to_string(Right.space.size()) + ").");
    }

    // Do the action
    return Left*Right;
}
}
}

#endif
