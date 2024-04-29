#ifndef _LIELAB_FUNCTIONS_ACTIONS_HPP
#define _LIELAB_FUNCTIONS_ACTIONS_HPP

#include "exp.hpp"

namespace lielab
{
namespace functions
{
lielab::domain::hmlie left_product(const lielab::domain::hmlie & Left, const lielab::domain::hmlie & Right)
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
    // TODO: Move this error checking to hmlie?
    if (Left.space.size() != Right.space.size())
    {
        throw SizeError("left_product: Spaces of hmlie must be the same size (" + std::to_string(Left.space.size()) + " != " + std::to_string(Right.space.size()) + ").");
    }

    // Do the action
    return Left*Right;
}
}
}

#endif
