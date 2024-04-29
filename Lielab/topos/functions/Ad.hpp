#ifndef _LIELAB_TOPOS_AD_HPP
#define _LIELAB_TOPOS_AD_HPP

#include "../../domain.hpp"
#include "../../functions.hpp"

namespace Lielab
{
namespace topos
{
Lielab::domain::CompositeAlgebra Ad(const Lielab::domain::CompositeManifold & A, const Lielab::domain::CompositeAlgebra & b)
{
    /*!
    * CompositeAlgebra Ad overload.
    */

    Lielab::domain::CompositeAlgebra out;

    for (int ii = 0; ii < b.space.size(); ii++)
    {
        const size_t ind = b.space[ii].index();
        if (ind == Lielab::domain::CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<Lielab::domain::CN>(A.space[ii]),
                                                        std::get<Lielab::domain::cn>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_gl)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<Lielab::domain::GL>(A.space[ii]),
                                                        std::get<Lielab::domain::gl>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<Lielab::domain::RN>(A.space[ii]),
                                                        std::get<Lielab::domain::rn>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<Lielab::domain::SE>(A.space[ii]),
                                                        std::get<Lielab::domain::se>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<Lielab::domain::SO>(A.space[ii]),
                                                        std::get<Lielab::domain::so>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<Lielab::domain::SP>(A.space[ii]),
                                                      std::get<Lielab::domain::sp>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<Lielab::domain::SU>(A.space[ii]),
                                                      std::get<Lielab::domain::su>(b.space[ii])));
        }
    }

    return out;
}
}
}

#endif
