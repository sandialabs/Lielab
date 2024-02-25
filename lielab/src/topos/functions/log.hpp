#ifndef _LIELAB_TOPOS_LOG_HPP
#define _LIELAB_TOPOS_LOG_HPP

namespace lielab
{
namespace topos
{
lielab::domain::halie log(const lielab::domain::hmlie &M)
{
    /*!
     * hmlie logarithm overload
     */

    lielab::domain::halie out;
    
    for (int ii = 0; ii < M.space.size(); ii++)
    {
        const size_t ind = M.space[ii].index();
        if (ind == lielab::domain::hmlie::INDEX_GL)
        {
            out.space.push_back(lielab::functions::log(std::get<lielab::domain::GL>(M.space[ii])));
        }
        else if (ind == lielab::domain::hmlie::INDEX_RN)
        {
            out.space.push_back(lielab::functions::log(std::get<lielab::domain::RN>(M.space[ii])));
        }
        else if (ind == lielab::domain::hmlie::INDEX_SO)
        {
            out.space.push_back(lielab::functions::log(std::get<lielab::domain::SO>(M.space[ii])));
        }
        else if (ind == lielab::domain::hmlie::INDEX_SP)
        {
            out.space.push_back(lielab::functions::log(std::get<lielab::domain::SP>(M.space[ii])));
        }
        else if (ind == lielab::domain::hmlie::INDEX_SU)
        {
            out.space.push_back(lielab::functions::log(std::get<lielab::domain::SU>(M.space[ii])));
        }
    }

    return out;
}
}
}

#endif
