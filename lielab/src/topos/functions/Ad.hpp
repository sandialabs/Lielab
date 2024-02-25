#ifndef _LIELAB_TOPOS_AD_HPP
#define _LIELAB_TOPOS_AD_HPP

namespace lielab
{
namespace topos
{
lielab::domain::halie Ad(const lielab::domain::hmlie & A, const lielab::domain::halie & b)
{
    /*!
    * halie Ad overload.
    */

    lielab::domain::halie out;

    for (int ii = 0; ii < b.space.size(); ii++)
    {
        const size_t ind = b.space[ii].index();
        if (ind == lielab::domain::halie::INDEX_gl)
        {
            out.space.push_back(lielab::functions::Ad(std::get<lielab::domain::GL>(A.space[ii]),
                                                        std::get<lielab::domain::gl>(b.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_rn)
        {
            out.space.push_back(lielab::functions::Ad(std::get<lielab::domain::RN>(A.space[ii]),
                                                        std::get<lielab::domain::rn>(b.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_so)
        {
            out.space.push_back(lielab::functions::Ad(std::get<lielab::domain::SO>(A.space[ii]),
                                                        std::get<lielab::domain::so>(b.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_sp)
        {
            out.space.push_back(lielab::functions::Ad(std::get<lielab::domain::SP>(A.space[ii]),
                                                      std::get<lielab::domain::sp>(b.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_su)
        {
            out.space.push_back(lielab::functions::Ad(std::get<lielab::domain::SU>(A.space[ii]),
                                                      std::get<lielab::domain::su>(b.space[ii])));
        }
    }

    return out;
}
}
}

#endif
