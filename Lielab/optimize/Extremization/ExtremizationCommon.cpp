#include "ExtremizationCommon.hpp"

#include <functional>

namespace Lielab::optimize
{

EuclideanExtremizationSystem::EuclideanExtremizationSystem(EuclideanExtremizationSystem_objective_t objective_)
{
    this->objective = objective_;
}

}
