#ifndef _LIELAB_TOPOS_DEXPINV_HPP
#define _LIELAB_TOPOS_DEXPINV_HPP

#include "../../domain.hpp"
#include "../../functions.hpp"

namespace Lielab
{
namespace topos
{
Lielab::domain::CompositeAlgebra dexpinv(const Lielab::domain::CompositeAlgebra & a, const Lielab::domain::CompositeAlgebra & b)
{
    /*!
    * CompositeAlgebra dexpinv overload
    */

    Lielab::domain::CompositeAlgebra out;

    for (int ii = 0; ii < a.space.size(); ii++)
    {
        const size_t ind = a.space[ii].index();
        if (ind == Lielab::domain::CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::cn>(a.space[ii]),
                                                           std::get<Lielab::domain::cn>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_gl)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::gl>(a.space[ii]),
                                                           std::get<Lielab::domain::gl>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::glc>(a.space[ii]),
                                                           std::get<Lielab::domain::glc>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::rn>(a.space[ii]),
                                                           std::get<Lielab::domain::rn>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::se>(a.space[ii]),
                                                           std::get<Lielab::domain::se>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::so>(a.space[ii]),
                                                           std::get<Lielab::domain::so>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::sp>(a.space[ii]),
                                                           std::get<Lielab::domain::sp>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(Lielab::functions::dexpinv(std::get<Lielab::domain::su>(a.space[ii]),
                                                           std::get<Lielab::domain::su>(b.space[ii])));
        }
    }

    return out;
}
}
}

#endif
