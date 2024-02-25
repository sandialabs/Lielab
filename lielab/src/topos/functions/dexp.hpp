#ifndef _LIELAB_TOPOS_DEXP_HPP
#define _LIELAB_TOPOS_DEXP_HPP

namespace lielab
{
namespace topos
{
lielab::domain::halie dexp(const lielab::domain::halie & a, const lielab::domain::halie & b, const size_t order = 5)
{
    /*!
    * halie d exponential overload.
    */

    lielab::domain::halie out;

    for (int ii = 0; ii < a.space.size(); ii++)
    {
        const size_t ind = a.space[ii].index();
        if (ind == lielab::domain::halie::INDEX_gl)
        {
            out.space.push_back(lielab::functions::dexp(std::get<lielab::domain::gl>(a.space[ii]),
                                                        std::get<lielab::domain::gl>(b.space[ii]), order));
        }
        else if (ind == lielab::domain::halie::INDEX_rn)
        {
            out.space.push_back(lielab::functions::dexp(std::get<lielab::domain::rn>(a.space[ii]),
                                                        std::get<lielab::domain::rn>(b.space[ii]), order));
        }
        else if (ind == lielab::domain::halie::INDEX_so)
        {
            out.space.push_back(lielab::functions::dexp(std::get<lielab::domain::so>(a.space[ii]),
                                                        std::get<lielab::domain::so>(b.space[ii]), order));
        }
        else if (ind == lielab::domain::halie::INDEX_sp)
        {
            out.space.push_back(lielab::functions::dexp(std::get<lielab::domain::sp>(a.space[ii]),
                                                        std::get<lielab::domain::sp>(b.space[ii]), order));
        }
        else if (ind == lielab::domain::halie::INDEX_su)
        {
            out.space.push_back(lielab::functions::dexp(std::get<lielab::domain::su>(a.space[ii]),
                                                        std::get<lielab::domain::su>(b.space[ii]), order));
        }
    }

    return out;
}
}
}

#endif
