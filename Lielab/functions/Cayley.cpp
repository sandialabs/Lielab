#include "Cayley.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

template <>
Lielab::domain::CN Cayley(const Lielab::domain::cn & a)
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
Lielab::domain::GLC Cayley(const Lielab::domain::glc & a)
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
Lielab::domain::SU Cayley(const Lielab::domain::su & a)
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

Lielab::domain::CompositeGroup Cayley(const Lielab::domain::CompositeAlgebra & la)
{
    /*!
    * CompositeAlgebra cayley overload.
    */

    Lielab::domain::CompositeGroup out;

    for (size_t ii = 0; ii < la.space.size(); ii++)
    {
        const size_t ind = la.space[ii].index();
        if (ind == Lielab::domain::CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<Lielab::domain::cn>(la.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<Lielab::domain::glr>(la.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<Lielab::domain::glc>(la.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<Lielab::domain::rn>(la.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<Lielab::domain::se>(la.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<Lielab::domain::so>(la.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<Lielab::domain::sp>(la.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(Lielab::functions::Cayley(std::get<Lielab::domain::su>(la.space[ii])));
        }
    }

    return out;
}

Lielab::domain::CompositeAlgebra dCayleyinv(const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b)
{
    /*!
    * CompositeAlgebra dCayleyinv overload
    */

    Lielab::domain::CompositeAlgebra out;

    for (size_t ii = 0; ii < a.space.size(); ii++)
    {
        const size_t ind = a.space[ii].index();
        if (ind == Lielab::domain::CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<Lielab::domain::cn>(a.space[ii]),
                                                              std::get<Lielab::domain::cn>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<Lielab::domain::glr>(a.space[ii]),
                                                              std::get<Lielab::domain::glr>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<Lielab::domain::glc>(a.space[ii]),
                                                              std::get<Lielab::domain::glc>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<Lielab::domain::rn>(a.space[ii]),
                                                              std::get<Lielab::domain::rn>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<Lielab::domain::se>(a.space[ii]),
                                                              std::get<Lielab::domain::se>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<Lielab::domain::so>(a.space[ii]),
                                                              std::get<Lielab::domain::so>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<Lielab::domain::sp>(a.space[ii]),
                                                              std::get<Lielab::domain::sp>(b.space[ii])));
        }
        else if (ind == Lielab::domain::CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(Lielab::functions::dCayleyinv(std::get<Lielab::domain::su>(a.space[ii]),
                                                              std::get<Lielab::domain::su>(b.space[ii])));
        }
    }

    return out;
}

}
