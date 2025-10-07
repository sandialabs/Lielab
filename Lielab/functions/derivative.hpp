#ifndef LIELAB_FUNCTIONS_DERIVATIVE_HPP
#define LIELAB_FUNCTIONS_DERIVATIVE_HPP

#include "Lielab/domain.hpp"

#include <functional>

namespace Lielab::functions
{

Lielab::domain::CompositeAlgebra forward_difference(std::function<Lielab::domain::CompositeGroup(const double)> & fun, const double t, const double dt);
Lielab::domain::CompositeAlgebra backward_difference(std::function<Lielab::domain::CompositeGroup(const double)> & fun, const double t, const double dt);
Lielab::domain::CompositeAlgebra central_difference(std::function<Lielab::domain::CompositeGroup(const double)> & fun, const double t, const double dt);

}

#endif
