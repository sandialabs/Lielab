#ifndef EIGEN_MPL2_ONLY
#define EIGEN_MPL2_ONLY
#endif
#define CATCH_CONFIG_MAIN

// Standard includes
#include <cmath>
#include <iostream>
#include <vector>

// Additional includes (extra code, etc)
#include <catch2/catch_all.hpp>
#include <Eigen/Core>
#include <lielab>

// Include test utils specific for lielab
#include "test_utils.hpp"

// Test includes
#include "test_domain.cpp"
#include "test_dynamics.cpp"
#include "test_functions.cpp"
#include "test_kinematics.cpp"
#include "test_optim.cpp"
// #include "test_lielab.cpp"
#include "test_topos.cpp"
#include "test_transform.cpp"
