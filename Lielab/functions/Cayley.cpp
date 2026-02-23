#include "Cayley.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/testing.hpp"

namespace Lielab::functions
{

template <>
Lielab::domain::CompositeGroup cay(const Lielab::domain::CompositeAlgebra& x)
{
    /*!
    * CompositeAlgebra cayley overload.
    */

    using namespace Lielab::domain;

    CompositeGroup out;
    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(cay(_element));
        }, element);
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra cayinv(const Lielab::domain::CompositeGroup& g)
{
    /*!
    * CompositeAlgebra cayley overload.
    */

    using namespace Lielab::domain;

    CompositeAlgebra out;
    for (const auto& element : g.point)
    {
        std::visit([&](const auto& _element)
        {
            out.point.push_back(cayinv(_element));
        }, element);
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dcay(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y)
{
    /*!
    * CompositeAlgebra dcay overload
    */

    using Lielab::domain::CompositeAlgebra;
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(x, y), "Unable to take dcay of topologically inconsistent algebras: (" + x.to_string() + ") !≅ (" + y.to_string() + ").");

    CompositeAlgebra out;
    int index = 0;
    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            out.point.push_back(dcay(_element, std::get<other_t>(y.point[index])));
        }, element);
        index++;
    }

    return out;
}

template <>
Lielab::domain::CompositeAlgebra dcayinv(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y)
{
    /*!
    * CompositeAlgebra dcayinv overload
    */

    using Lielab::domain::CompositeAlgebra;
    using Lielab::testing::check_topology;

    lielab_assert(check_topology(x, y), "Unable to take dcayinv of topologically inconsistent algebras: (" + x.to_string() + ") !≅ (" + y.to_string() + ").");

    CompositeAlgebra out;
    int index = 0;
    for (const auto& element : x.point)
    {
        std::visit([&](const auto& _element)
        {
            using other_t = std::decay_t<decltype(_element)>;
            out.point.push_back(dcayinv(_element, std::get<other_t>(y.point[index])));
        }, element);
        index++;
    }

    return out;
}

}
