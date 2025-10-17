#include "Ad.hpp"

#include "littlead.hpp"
#include "exp.hpp"

#include "Lielab/domain.hpp"

#include <cassert>

namespace Lielab::functions
{

Lielab::domain::CompositeAlgebra Ad(const Lielab::domain::CompositeGroup& A, const Lielab::domain::CompositeAlgebra& b)
{
    /*!
    * CompositeAlgebra Ad overload.
    */

    using namespace Lielab::domain;

    CompositeAlgebra out;

    for (size_t ii = 0; ii < b.space.size(); ii++)
    {
        const size_t ind = b.space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<CN>(A.space[ii]),
                                                        std::get<cn>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<GLR>(A.space[ii]),
                                                        std::get<glr>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<GLC>(A.space[ii]),
                                                        std::get<glc>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<RN>(A.space[ii]),
                                                        std::get<rn>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<SE>(A.space[ii]),
                                                        std::get<se>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<SO>(A.space[ii]),
                                                        std::get<so>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<SP>(A.space[ii]),
                                                      std::get<sp>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(Lielab::functions::Ad(std::get<SU>(A.space[ii]),
                                                      std::get<su>(b.space[ii])));
        }
    }

    return out;
}

}
