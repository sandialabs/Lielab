#ifndef _LIELAB_FUNCTIONS_ACTIONS_HPP
#define _LIELAB_FUNCTIONS_ACTIONS_HPP

#include "exp.hpp"

namespace lielab
{
namespace functions
{
lielab::domain::hmlie left_exp_default(const lielab::domain::halie & left, const lielab::domain::hmlie & Right)
{
    /*!
     * Default exponential left action of an algebra on a manifold
     *
     * @param[in] left
     * @param[in] Right
     * @param[out] out
     * 
     * Author / Date: Sparapany / 2022
     */

    // Simple error checking on the inputs
    if (left.space.size() != Right.space.size())
    {
        throw SizeError("left_exp_default: Spaces of halie must be the same size (" + std::to_string(left.space.size()) + " != " + std::to_string(Right.space.size()) + ").");
    }
    
    // Begin by taking the exponential of the algebra
    lielab::domain::hmlie Left;

    for (int ii = 0; ii < left.space.size(); ii++)
    {
        const size_t ind = left.space[ii].index();
        if (ind == lielab::domain::halie::INDEX_gl)
        {
            Left.space.push_back(lielab::functions::exp(std::get<lielab::domain::gl>(left.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_rn)
        {
            Left.space.push_back(lielab::functions::exp(std::get<lielab::domain::rn>(left.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_so)
        {
            Left.space.push_back(lielab::functions::exp(std::get<lielab::domain::so>(left.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_sp)
        {
            Left.space.push_back(lielab::functions::exp(std::get<lielab::domain::sp>(left.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_su)
        {
            Left.space.push_back(lielab::functions::exp(std::get<lielab::domain::su>(left.space[ii])));
        }
    }

    // After error checking, the resulting Left should have the same size
    // as the input right, otherwise something is wrong in the loop above.
    assert(Left.space.size() == Right.space.size());

    // Do the action on Right
    lielab::domain::hmlie out;

    for (int ii = 0; ii < Left.space.size(); ii++)
    {
        const size_t ind1 = Left.space[ii].index();
        const size_t ind2 = Right.space[ii].index();
        if (ind1 == lielab::domain::hmlie::INDEX_GL)
        {
            out.space.push_back(std::get<lielab::domain::GL>(Left.space[ii]) * std::get<lielab::domain::GL>(Right.space[ii]));
        }
        else if (ind1 == lielab::domain::hmlie::INDEX_RN)
        {
            out.space.push_back(std::get<lielab::domain::RN>(Left.space[ii]) * std::get<lielab::domain::RN>(Right.space[ii]));
        }
        else if (ind1 == lielab::domain::hmlie::INDEX_SO)
        {
            out.space.push_back(std::get<lielab::domain::SO>(Left.space[ii]) * std::get<lielab::domain::SO>(Right.space[ii]));
        }
        else if (ind1 == lielab::domain::hmlie::INDEX_SP)
        {
            out.space.push_back(std::get<lielab::domain::SP>(Left.space[ii]) * std::get<lielab::domain::SP>(Right.space[ii]));
        }
        else if (ind1 == lielab::domain::hmlie::INDEX_SU)
        {
            out.space.push_back(std::get<lielab::domain::SU>(Left.space[ii]) * std::get<lielab::domain::SU>(Right.space[ii]));
        }
    }
    return out;
}
}
}

#endif
