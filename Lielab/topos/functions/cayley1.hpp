#ifndef _LIELAB_TOPOS_CAYLEY1_HPP
#define _LIELAB_TOPOS_CAYLEY1_HPP

#include "../../domain.hpp"
#include "../../functions.hpp"

namespace Lielab
{
    namespace topos
    {
        Lielab::domain::CompositeManifold cayley1(const Lielab::domain::CompositeAlgebra & la)
        {
            /*!
            * CompositeAlgebra cayley overload.
            */

            Lielab::domain::CompositeManifold out;

            for (int ii = 0; ii < la.space.size(); ii++)
            {
                const size_t ind = la.space[ii].index();
                if (ind == Lielab::domain::CompositeAlgebra::INDEX_cn)
                {
                    out.space.push_back(Lielab::functions::cayley1(std::get<Lielab::domain::cn>(la.space[ii])));
                }
                else if (ind == Lielab::domain::CompositeAlgebra::INDEX_gl)
                {
                    out.space.push_back(Lielab::functions::cayley1(std::get<Lielab::domain::gl>(la.space[ii])));
                }
                else if (ind == Lielab::domain::CompositeAlgebra::INDEX_rn)
                {
                    out.space.push_back(Lielab::functions::cayley1(std::get<Lielab::domain::rn>(la.space[ii])));
                }
                else if (ind == Lielab::domain::CompositeAlgebra::INDEX_se)
                {
                    out.space.push_back(Lielab::functions::cayley1(std::get<Lielab::domain::se>(la.space[ii])));
                }
                else if (ind == Lielab::domain::CompositeAlgebra::INDEX_so)
                {
                    out.space.push_back(Lielab::functions::cayley1(std::get<Lielab::domain::so>(la.space[ii])));
                }
                else if (ind == Lielab::domain::CompositeAlgebra::INDEX_sp)
                {
                    out.space.push_back(Lielab::functions::cayley1(std::get<Lielab::domain::sp>(la.space[ii])));
                }
                else if (ind == Lielab::domain::CompositeAlgebra::INDEX_su)
                {
                    out.space.push_back(Lielab::functions::cayley1(std::get<Lielab::domain::su>(la.space[ii])));
                }
            }

            return out;
        }
    }
}

#endif
