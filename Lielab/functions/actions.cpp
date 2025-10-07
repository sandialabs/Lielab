#include "actions.hpp"

#include "exp.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/utils/Error.hpp"

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
    if (g.space.size() != y.space.size())
    {
        throw Lielab::utils::Error("left_Lie_group_action: CompositeGroup and CompositeManifold must be the same size (" + std::to_string(g.space.size()) + " != " + std::to_string(y.space.size()) + ").");
    }

    // Do the action
    CompositeManifold out;

    for (size_t ii = 0; ii < g.space.size(); ii++)
    {
        const size_t indg = g.space[ii].index();
        const size_t indy = y.space[ii].index();
        if (indg == CompositeGroup::INDEX_CN && indy == CompositeManifold::INDEX_CN)
        {
            out.space.push_back(std::get<Lielab::domain::CN>(g.space[ii]) * std::get<Lielab::domain::CN>(y.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_GLR && indy == CompositeManifold::INDEX_GLR)
        {
            out.space.push_back(std::get<Lielab::domain::GLR>(g.space[ii]) * std::get<Lielab::domain::GLR>(y.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_GLC && indy == CompositeManifold::INDEX_GLC)
        {
            out.space.push_back(std::get<Lielab::domain::GLC>(g.space[ii]) * std::get<Lielab::domain::GLC>(y.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_RN && indy == CompositeManifold::INDEX_RN)
        {
            out.space.push_back(std::get<Lielab::domain::RN>(g.space[ii]) * std::get<Lielab::domain::RN>(y.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SE && indy == CompositeManifold::INDEX_SE)
        {
            out.space.push_back(std::get<Lielab::domain::SE>(g.space[ii]) * std::get<Lielab::domain::SE>(y.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SO && indy == CompositeManifold::INDEX_SO)
        {
            out.space.push_back(std::get<Lielab::domain::SO>(g.space[ii]) * std::get<Lielab::domain::SO>(y.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SP && indy == CompositeManifold::INDEX_SP)
        {
            out.space.push_back(std::get<Lielab::domain::SP>(g.space[ii]) * std::get<Lielab::domain::SP>(y.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SU && indy == CompositeManifold::INDEX_SU)
        {
            out.space.push_back(std::get<Lielab::domain::SU>(g.space[ii]) * std::get<Lielab::domain::SU>(y.space[ii]));
        }
        else
        {
            throw Lielab::utils::Error("left_Lie_group_action: Unknown action with given structure.");
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
    if (g.space.size() != y.space.size())
    {
        throw Lielab::utils::Error("right_Lie_group_action: CompositeGroup and CompositeManifold must be the same size (" + std::to_string(g.space.size()) + " != " + std::to_string(y.space.size()) + ").");
    }

    // Do the action
    CompositeManifold out;

    for (size_t ii = 0; ii < g.space.size(); ii++)
    {
        const size_t indg = g.space[ii].index();
        const size_t indy = y.space[ii].index();
        if (indg == CompositeGroup::INDEX_CN && indy == CompositeManifold::INDEX_CN)
        {
            out.space.push_back(std::get<Lielab::domain::CN>(y.space[ii]) * std::get<Lielab::domain::CN>(g.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_GLR && indy == CompositeManifold::INDEX_GLR)
        {
            out.space.push_back(std::get<Lielab::domain::GLR>(y.space[ii]) * std::get<Lielab::domain::GLR>(g.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_GLC && indy == CompositeManifold::INDEX_GLC)
        {
            out.space.push_back(std::get<Lielab::domain::GLC>(y.space[ii]) * std::get<Lielab::domain::GLC>(g.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_RN && indy == CompositeManifold::INDEX_RN)
        {
            out.space.push_back(std::get<Lielab::domain::RN>(y.space[ii]) * std::get<Lielab::domain::RN>(g.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SE && indy == CompositeManifold::INDEX_SE)
        {
            out.space.push_back(std::get<Lielab::domain::SE>(y.space[ii]) * std::get<Lielab::domain::SE>(g.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SO && indy == CompositeManifold::INDEX_SO)
        {
            out.space.push_back(std::get<Lielab::domain::SO>(y.space[ii]) * std::get<Lielab::domain::SO>(g.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SP && indy == CompositeManifold::INDEX_SP)
        {
            out.space.push_back(std::get<Lielab::domain::SP>(y.space[ii]) * std::get<Lielab::domain::SP>(g.space[ii]));
        }
        else if (indg == CompositeGroup::INDEX_SU && indy == CompositeManifold::INDEX_SU)
        {
            out.space.push_back(std::get<Lielab::domain::SU>(y.space[ii]) * std::get<Lielab::domain::SU>(g.space[ii]));
        }
        else
        {
            throw Lielab::utils::Error("right_Lie_group_action: Unknown action with given structure.");
        }
    }

    return out;
}

Lielab::domain::CompositeManifold left_Lie_algebra_action(const Lielab::domain::CompositeAlgebra& xi, const Lielab::domain::CompositeManifold& y)
{
    /*!
     * Default action by left product.
     *
     * Designed to just be a placeholder.
     *
     * Note: There is no tangent object so manifold outputs should be interacted with by get_matrix() or serialize().
     *
     * @param[in] xi
     * @param[in] y
     * @param[out] out g*y
     * 
     * Author / Date: Sparapany / 2025
     */

    using namespace Lielab::domain;

    // Simple error checking on the inputs
    if (xi.space.size() != y.space.size())
    {
        throw Lielab::utils::Error("left_Lie_algebra_action: CompositeAlgebra and CompositeManifold must be the same size (" + std::to_string(xi.space.size()) + " != " + std::to_string(y.space.size()) + ").");
    }

    // Do the daction
    CompositeManifold out;

    for (size_t ii = 0; ii < xi.space.size(); ii++)
    {
        const size_t indxi = xi.space[ii].index();
        const size_t indy = y.space[ii].index();

        if (indxi == CompositeAlgebra::INDEX_cn && indy == CompositeManifold::INDEX_CN)
        {
            const cn a = std::get<cn>(xi.space[ii]);
            const CN b = std::get<CN>(y.space[ii]);
            const CN dx = CN::from_complex_vector(a.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_glr && indy == CompositeManifold::INDEX_GLR)
        {
            const glr a = std::get<glr>(xi.space[ii]);
            const GLR b = std::get<GLR>(y.space[ii]);
            const GLR dx = GLR(a.data*b.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_glc && indy == CompositeManifold::INDEX_GLC)
        {
            const glc a = std::get<glc>(xi.space[ii]);
            const GLC b = std::get<GLC>(y.space[ii]);
            const GLC dx = GLC(a.data*b.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_rn && indy == CompositeManifold::INDEX_RN)
        {
            const rn a = std::get<rn>(xi.space[ii]);
            const RN b = std::get<RN>(y.space[ii]);
            const RN dx = RN::from_vector(a.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_se && indy == CompositeManifold::INDEX_SE)
        {
            const se a = std::get<se>(xi.space[ii]);
            const SE b = std::get<SE>(y.space[ii]);
            const SE dx = SE(a.data*b.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_so && indy == CompositeManifold::INDEX_SO)
        {
            const so a = std::get<so>(xi.space[ii]);
            const SO b = std::get<SO>(y.space[ii]);
            const SO dx = SO(a.data*b.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_sp && indy == CompositeManifold::INDEX_SP)
        {
            const sp a = std::get<sp>(xi.space[ii]);
            const SP b = std::get<SP>(y.space[ii]);
            const SP dx = SP(a.data*b.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_su && indy == CompositeManifold::INDEX_SU)
        {
            const su a = std::get<su>(xi.space[ii]);
            const SU b = std::get<SU>(y.space[ii]);
            const SU dx = SU(a.data*b.data);
            out.space.push_back(dx);
        }
        else
        {
            throw Lielab::utils::Error("left_Lie_algebra_action: Unknown action with given structure.");
        }
    }

    return out;
}

Lielab::domain::CompositeManifold right_Lie_algebra_action(const Lielab::domain::CompositeAlgebra& xi, const Lielab::domain::CompositeManifold& y)
{
    /*!
     * Default action by right product.
     *
     * Designed to just be a placeholder.
     *
     * Note: There is no tangent object so manifold outputs should be interacted with by get_matrix() or serialize().
     *
     * @param[in] xi
     * @param[in] y
     * @param[out] out g*y
     * 
     * Author / Date: Sparapany / 2025
     */

    using namespace Lielab::domain;

    // Simple error checking on the inputs
    if (xi.space.size() != y.space.size())
    {
        throw Lielab::utils::Error("right_Lie_algebra_action: CompositeAlgebra and CompositeManifold must be the same size (" + std::to_string(xi.space.size()) + " != " + std::to_string(y.space.size()) + ").");
    }

    // Do the daction
    CompositeManifold out;

    for (size_t ii = 0; ii < xi.space.size(); ii++)
    {
        const size_t indxi = xi.space[ii].index();
        const size_t indy = y.space[ii].index();

        if (indxi == CompositeAlgebra::INDEX_cn && indy == CompositeManifold::INDEX_CN)
        {
            const cn a = std::get<cn>(xi.space[ii]);
            const CN b = std::get<CN>(y.space[ii]);
            const CN dx = CN::from_complex_vector(a.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_glr && indy == CompositeManifold::INDEX_GLR)
        {
            const glr a = std::get<glr>(xi.space[ii]);
            const GLR b = std::get<GLR>(y.space[ii]);
            const GLR dx = GLR(b.data*a.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_glc && indy == CompositeManifold::INDEX_GLC)
        {
            const glc a = std::get<glc>(xi.space[ii]);
            const GLC b = std::get<GLC>(y.space[ii]);
            const GLC dx = GLC(b.data*a.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_rn && indy == CompositeManifold::INDEX_RN)
        {
            const rn a = std::get<rn>(xi.space[ii]);
            const RN b = std::get<RN>(y.space[ii]);
            const RN dx = RN::from_vector(a.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_se && indy == CompositeManifold::INDEX_SE)
        {
            const se a = std::get<se>(xi.space[ii]);
            const SE b = std::get<SE>(y.space[ii]);
            const SE dx = SE(b.data*a.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_so && indy == CompositeManifold::INDEX_SO)
        {
            const so a = std::get<so>(xi.space[ii]);
            const SO b = std::get<SO>(y.space[ii]);
            const SO dx = SO(b.data*a.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_sp && indy == CompositeManifold::INDEX_SP)
        {
            const sp a = std::get<sp>(xi.space[ii]);
            const SP b = std::get<SP>(y.space[ii]);
            const SP dx = SP(b.data*a.data);
            out.space.push_back(dx);
        }
        else if (indxi == CompositeAlgebra::INDEX_su && indy == CompositeManifold::INDEX_SU)
        {
            const su a = std::get<su>(xi.space[ii]);
            const SU b = std::get<SU>(y.space[ii]);
            const SU dx = SU(b.data*a.data);
            out.space.push_back(dx);
        }
        else
        {
            throw Lielab::utils::Error("right_Lie_algebra_action: Unknown action with given structure.");
        }
    }

    return out;
}

}
