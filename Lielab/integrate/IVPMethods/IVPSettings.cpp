#include "IVPSettings.hpp"

namespace Lielab::integrate
{

EuclideanIVPSystem::EuclideanIVPSystem(EuclideanIVP_vectorfield_t vf)
{
    this->vectorfield = vf;
}

HomogeneousIVPSystem::HomogeneousIVPSystem(HomogeneousIVP_vectorfield_t vf)
{
    this->vectorfield = vf;
}

}
