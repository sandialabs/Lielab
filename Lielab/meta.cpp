#include <Eigen/Core>

#include "meta.hpp"

#include <string>

namespace Lielab
{

std::string get_simd_info()
{
    return Eigen::SimdInstructionSetsInUse();
}

std::string get_eigen_info()
{
    return "Eigen3 Version: " + std::to_string(EIGEN_MAJOR_VERSION) + "." + std::to_string(EIGEN_MINOR_VERSION) + "." + std::to_string(EIGEN_PATCH_VERSION);
}

}
