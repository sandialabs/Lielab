#ifndef LIELAB_FUNCTIONS_ACTIONS_HPP
#define LIELAB_FUNCTIONS_ACTIONS_HPP

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

Lielab::domain::CompositeManifold left_Lie_group_action(const Lielab::domain::CompositeGroup& g, const Lielab::domain::CompositeManifold& y);
Lielab::domain::CompositeManifold right_Lie_group_action(const Lielab::domain::CompositeGroup& g, const Lielab::domain::CompositeManifold& y);

}

#endif
