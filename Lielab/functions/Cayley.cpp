#include "Cayley.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <>
Lielab::domain::CN Cayley(const Lielab::domain::cn& a)
{
    /*
    * Cayley overload for cn
    *
    * Needed since cn is complex
    */

    const size_t shape = a.get_shape();

    const Eigen::MatrixXcd m = a.get_matrix();
    const Eigen::MatrixXcd Id = Eigen::MatrixXcd::Identity(shape, shape);

    return (Id + m/2.0)*(Id - m/2.0).inverse();
}

template <>
Lielab::domain::GLC Cayley(const Lielab::domain::glc& a)
{
    /*
    * Cayley overload for glc
    *
    * Needed since glc is complex
    */

    const size_t shape = a.get_shape();

    const Eigen::MatrixXcd m = a.get_matrix();
    const Eigen::MatrixXcd Id = Eigen::MatrixXcd::Identity(shape, shape);

    return (Id + m/2.0)*(Id - m/2.0).inverse();
}

template <>
Lielab::domain::SU Cayley(const Lielab::domain::su& a)
{
    /*
    * Cayley overload for su
    *
    * Needed since su is complex
    */

    const size_t shape = a.get_shape();

    const Eigen::MatrixXcd m = a.get_matrix();
    const Eigen::MatrixXcd Id = Eigen::MatrixXcd::Identity(shape, shape);

    return (Id + m/2.0)*(Id - m/2.0).inverse();
}

Lielab::domain::CompositeGroup Cayley(const Lielab::domain::CompositeAlgebra& la)
{
    /*!
    * CompositeAlgebra cayley overload.
    */

    using namespace Lielab::domain;

    CompositeGroup out;

    for (size_t ii = 0; ii < la.space.size(); ii++)
    {
        const size_t ind = la.space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<cn>(la.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<glr>(la.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<glc>(la.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<rn>(la.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<se>(la.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<so>(la.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<sp>(la.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<su>(la.space[ii])));
        }
    }

    return out;
}

Lielab::domain::CompositeAlgebra dCayleyinv(const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b)
{
    /*!
    * CompositeAlgebra dCayleyinv overload
    */

    using namespace Lielab::domain;

    CompositeAlgebra out;

    for (size_t ii = 0; ii < a.space.size(); ii++)
    {
        const size_t ind = a.space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<cn>(a.space[ii]),
                                                              std::get<cn>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<glr>(a.space[ii]),
                                                              std::get<glr>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<glc>(a.space[ii]),
                                                              std::get<glc>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<rn>(a.space[ii]),
                                                              std::get<rn>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<se>(a.space[ii]),
                                                              std::get<se>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<so>(a.space[ii]),
                                                              std::get<so>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<sp>(a.space[ii]),
                                                              std::get<sp>(b.space[ii])));
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<su>(a.space[ii]),
                                                              std::get<su>(b.space[ii])));
        }
    }

    return out;
}

}
