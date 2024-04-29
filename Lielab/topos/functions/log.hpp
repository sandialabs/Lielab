#ifndef _LIELAB_TOPOS_LOG_HPP
#define _LIELAB_TOPOS_LOG_HPP

#include "../../domain.hpp"
#include "../../functions.hpp"

namespace Lielab
{
namespace topos
{

Lielab::domain::CompositeAlgebra log(const Lielab::domain::CompositeManifold &M)
{
    /*!
     * CompositeManifold logarithm overload
     */

    Lielab::domain::CompositeAlgebra out;
    
    for (int ii = 0; ii < M.space.size(); ii++)
    {
        const size_t ind = M.space[ii].index();
        if (ind == Lielab::domain::CompositeManifold::INDEX_CN)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::CN>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeManifold::INDEX_GL)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::GL>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeManifold::INDEX_GLC)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::GLC>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeManifold::INDEX_RN)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::RN>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeManifold::INDEX_SE)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::SE>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeManifold::INDEX_SO)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::SO>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeManifold::INDEX_SP)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::SP>(M.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeManifold::INDEX_SU)
        {
            out.space.push_back(Lielab::functions::log(std::get<Lielab::domain::SU>(M.space[ii])));
        }
    }

    return out;
}

}
}

#endif
