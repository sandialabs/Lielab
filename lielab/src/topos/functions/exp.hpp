#ifndef _LIELAB_TOPOS_EXP_HPP
#define _LIELAB_TOPOS_EXP_HPP

namespace lielab
{
namespace topos
{

lielab::domain::hmlie exp(const lielab::domain::halie & la)
{
    /*!
    * halie exponential overload.
    */

    lielab::domain::hmlie out;

    for (int ii = 0; ii < la.space.size(); ii++)
    {
        const size_t ind = la.space[ii].index();
        if (ind == lielab::domain::halie::INDEX_gl)
        {
            out.space.push_back(lielab::functions::exp(std::get<lielab::domain::gl>(la.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_rn)
        {
            out.space.push_back(lielab::functions::exp(std::get<lielab::domain::rn>(la.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_se)
        {
            out.space.push_back(lielab::functions::exp(std::get<lielab::domain::se>(la.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_so)
        {
            out.space.push_back(lielab::functions::exp(std::get<lielab::domain::so>(la.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_sp)
        {
            out.space.push_back(lielab::functions::exp(std::get<lielab::domain::sp>(la.space[ii])));
        }
        else if (ind == lielab::domain::halie::INDEX_su)
        {
            out.space.push_back(lielab::functions::exp(std::get<lielab::domain::su>(la.space[ii])));
        }
    }

    return out;
}

}
}

#endif
