#ifndef _LIELAB_DYNAMICS_NEWTON_HPP
#define _LIELAB_DYNAMICS_NEWTON_HPP

#include "../../domain"

namespace lielab
{
namespace dynamics
{

lielab::domain::rn newton_eom(const Eigen::VectorXd force, const double mass)
{
    /*!
     * Equations-of-motion on Newton's law.
     *
     * a = m/f
     * 
     * @param[in] force
     * @param[in] mass
     * @param[out] out
     */

    return lielab::domain::rn{force/mass};
}

lielab::domain::RN newton_action(const lielab::domain::rn g, const lielab::domain::RN M)
{
    /*!
     * Action of Newton's law.
     *
     * Addition.
     *
     * @param[in] g
     * @param[in] M
     * @param[out] M2
     */

    return lielab::domain::RN(g._data + M._data);
}

}
}

#endif
