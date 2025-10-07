#ifndef LIELAB_FUNCTIONS_ACTIONS_HPP
#define LIELAB_FUNCTIONS_ACTIONS_HPP

#include "exp.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

Lielab::domain::CompositeManifold left_Lie_group_action(const Lielab::domain::CompositeGroup& g, const Lielab::domain::CompositeManifold& y);
Lielab::domain::CompositeManifold right_Lie_group_action(const Lielab::domain::CompositeGroup& g, const Lielab::domain::CompositeManifold& y);

Lielab::domain::CompositeManifold left_Lie_algebra_action(const Lielab::domain::CompositeAlgebra& xi, const Lielab::domain::CompositeManifold& y);
Lielab::domain::CompositeManifold right_Lie_algebra_action(const Lielab::domain::CompositeAlgebra& xi, const Lielab::domain::CompositeManifold& y);

}

#endif
