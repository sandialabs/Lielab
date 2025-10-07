#include "Ad.hpp"

#include "littlead.hpp"
#include "exp.hpp"

#include "Lielab/domain.hpp"

#include <cassert>

namespace Lielab::functions
{

Lielab::domain::CompositeAlgebra Ad(const Lielab::domain::CompositeGroup & A, const Lielab::domain::CompositeAlgebra & b)
{
    /*!
    * CompositeAlgebra Ad overload.
    */

    Lielab::domain::CompositeAlgebra out;

    for (size_t ii = 0; ii < b.space.size(); ii++)
    {
        const size_t ind = b.space[ii].index();
        if (ind == Lielab::domain::CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<Lielab::domain::CN>(A.space[ii]),
                                                        std::get<Lielab::domain::cn>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<Lielab::domain::GLR>(A.space[ii]),
                                                        std::get<Lielab::domain::glr>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<Lielab::domain::GLC>(A.space[ii]),
                                                        std::get<Lielab::domain::glc>(b.space[ii])));
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
