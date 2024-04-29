#ifndef _LIELAB_TOPOS_DEXP_HPP
#define _LIELAB_TOPOS_DEXP_HPP

#include "../../domain.hpp"
#include "../../functions.hpp"

namespace Lielab
{
namespace topos
{
Lielab::domain::CompositeAlgebra dexp(const Lielab::domain::CompositeAlgebra & a, const Lielab::domain::CompositeAlgebra & b, const size_t order = 5)
{
    /*!
    * CompositeAlgebra d exponential overload.
    */

    Lielab::domain::CompositeAlgebra out;

    for (int ii = 0; ii < a.space.size(); ii++)
    {
        const size_t ind = a.space[ii].index();
        if (ind == Lielab::domain::CompositeAlgebra::INDEX_gl)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::gl>(a.space[ii]),
                                                        std::get<Lielab::domain::gl>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::rn>(a.space[ii]),
                                                        std::get<Lielab::domain::rn>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::se>(a.space[ii]),
                                                        std::get<Lielab::domain::se>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::so>(a.space[ii]),
                                                        std::get<Lielab::domain::so>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::sp>(a.space[ii]),
                                                        std::get<Lielab::domain::sp>(b.space[ii]), order));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(Lielab::functions::dexp(std::get<Lielab::domain::su>(a.space[ii]),
                                                        std::get<Lielab::domain::su>(b.space[ii]), order));
        }
    }

    return out;
}
}
}

#endif
