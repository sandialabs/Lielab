#include "commutator.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <>
Lielab::domain::CompositeAlgebra commutator(const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b)
{
    /*!
    * CompositeAlgebra commutator overload.
    */

    using namespace Lielab::domain;

    CompositeAlgebra out;

    for (size_t ii = 0; ii < a.space.size(); ii++)
    {
        const size_t ind = a.space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(commutator(std::get<cn>(a.space[ii]),
                                                        std::get<cn>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(commutator(std::get<glr>(a.space[ii]),
                                                        std::get<glr>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(commutator(std::get<glc>(a.space[ii]),
                                                        std::get<glc>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(commutator(std::get<rn>(a.space[ii]),
                                                        std::get<rn>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(commutator(std::get<se>(a.space[ii]),
                                                        std::get<se>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(commutator(std::get<so>(a.space[ii]),
                                                        std::get<so>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(commutator(std::get<sp>(a.space[ii]),
                                                        std::get<sp>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(commutator(std::get<su>(a.space[ii]),
                                                        std::get<su>(b.space[ii])));
        }
    }

    return out;
}

}
