#include "actions.hpp"

#include "exponential.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/testing/assertions.hpp"

#include <format>

namespace Lielab::functions
{

Lielab::domain::CompositeManifold left_Lie_group_action(const Lielab::domain::CompositeGroup& g, const Lielab::domain::CompositeManifold& y)
{
    /*!
     * Default action by left product.
     *
     * Designed to just be a placeholder.
     *
     * @param[in] g
     * @param[in] y
     * @param[out] out g*y
     * 
     * Author / Date: Sparapany / 2022
     */

    using namespace Lielab::domain;

    // Simple error checking on the inputs
    lielab_assert(g.point.size() == y.point.size(), std::format("CompositeGroup and CompositeManifold must be the same size ({} != {}).", g.point.size(), y.point.size()));

    // TODO: Somehow check the topology here

    // Do the action
    CompositeManifold out;

    for (size_t ii = 0; ii < g.point.size(); ii++)
    {
        const size_t indg = g.point[ii].index();
        const size_t indy = y.point[ii].index();
        if (indg == CompositeGroup::INDEX_CN && indy == CompositeManifold::INDEX_CN)
        {
            out.point.push_back(std::get<CN>(g.point[ii]) * std::get<CN>(y.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_GLC && indy == CompositeManifold::INDEX_GLC)
        {
            out.point.push_back(std::get<GLC>(g.point[ii]) * std::get<GLC>(y.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_GLR && indy == CompositeManifold::INDEX_GLR)
        {
            out.point.push_back(std::get<GLR>(g.point[ii]) * std::get<GLR>(y.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_RN && indy == CompositeManifold::INDEX_RN)
        {
            out.point.push_back(std::get<RN>(g.point[ii]) * std::get<RN>(y.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SE && indy == CompositeManifold::INDEX_SE)
        {
            out.point.push_back(std::get<SE>(g.point[ii]) * std::get<SE>(y.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SO && indy == CompositeManifold::INDEX_SO)
        {
            out.point.push_back(std::get<SO>(g.point[ii]) * std::get<SO>(y.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SP && indy == CompositeManifold::INDEX_SP)
        {
            out.point.push_back(std::get<SP>(g.point[ii]) * std::get<SP>(y.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SU && indy == CompositeManifold::INDEX_SU)
        {
            out.point.push_back(std::get<SU>(g.point[ii]) * std::get<SU>(y.point[ii]));
        }
        else
        {
            lielab_assert(false, "left_Lie_group_action: Unknown action with given structure.");
        }
    }

    return out;
}

Lielab::domain::CompositeManifold right_Lie_group_action(const Lielab::domain::CompositeGroup& g, const Lielab::domain::CompositeManifold& y)
{
    /*!
     * Action by right product.
     *
     * Designed to just be a placeholder.
     *
     * @param[in] g
     * @param[in] y
     * @param[out] out g*y
     * 
     * Author / Date: Sparapany / 2022
     */

    using namespace Lielab::domain;

    // Simple error checking on the inputs
    lielab_assert(g.point.size() == y.point.size(), std::format("CompositeGroup and CompositeManifold must be the same size ({} != {}).", g.point.size(), y.point.size()));

    // TODO: Somehow check the topology here

    // Do the action
    CompositeManifold out;

    for (size_t ii = 0; ii < g.point.size(); ii++)
    {
        const size_t indg = g.point[ii].index();
        const size_t indy = y.point[ii].index();
        if (indg == CompositeGroup::INDEX_CN && indy == CompositeManifold::INDEX_CN)
        {
            out.point.push_back(std::get<CN>(y.point[ii]) * std::get<CN>(g.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_GLC && indy == CompositeManifold::INDEX_GLC)
        {
            out.point.push_back(std::get<GLC>(y.point[ii]) * std::get<GLC>(g.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_GLR && indy == CompositeManifold::INDEX_GLR)
        {
            out.point.push_back(std::get<GLR>(y.point[ii]) * std::get<GLR>(g.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_RN && indy == CompositeManifold::INDEX_RN)
        {
            out.point.push_back(std::get<RN>(y.point[ii]) * std::get<RN>(g.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SE && indy == CompositeManifold::INDEX_SE)
        {
            out.point.push_back(std::get<SE>(y.point[ii]) * std::get<SE>(g.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SO && indy == CompositeManifold::INDEX_SO)
        {
            out.point.push_back(std::get<SO>(y.point[ii]) * std::get<SO>(g.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SP && indy == CompositeManifold::INDEX_SP)
        {
            out.point.push_back(std::get<SP>(y.point[ii]) * std::get<SP>(g.point[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SU && indy == CompositeManifold::INDEX_SU)
        {
            out.point.push_back(std::get<SU>(y.point[ii]) * std::get<SU>(g.point[ii]));
        }
        else
        {
            lielab_assert(false, "right_Lie_group_action: Unknown action with given structure.");
        }
    }

    return out;
}

}
