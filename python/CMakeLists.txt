project(pyLielab LANGUAGES CXX)

# Dependencies
if (NOT DEFINED Eigen3_FOUND)
    find_package(Eigen3 REQUIRED)
endif()

if (NOT DEFINED pybind11_FOUND)
    find_package(pybind11 REQUIRED)
endif()

if (NOT DEFINED Lielab_FOUND)
    find_package(Lielab REQUIRED)
endif()

# Add the pybind11 cpp file
SET( APP_PYD cppLielab )

# Root imports
file(GLOB pybind_LielabSubmodulesCPP "lielab/*.cpp")

# Domain submodule
file(GLOB pybind_LielabDomainCPP "lielab/domain/*.cpp")
file(GLOB pybind_LielabDomainLieAlgebrasCPP "lielab/domain/liealgebras/*.cpp")
file(GLOB pybind_LielabDomainLieGroupsCPP "lielab/domain/liegroups/*.cpp")
file(GLOB pybind_LielabDomainSmoothManifoldsCPP "lielab/domain/smoothmanifolds/*.cpp")

# Functions submodule
file(GLOB pybind_LielabFunctionsCPP "lielab/functions/*.cpp")

# Integrate submodule
file(GLOB pybind_LielabIntegrateCPP "lielab/integrate/*.cpp")
file(GLOB pybind_LielabIntegrateIVPMethodsCPP "lielab/integrate/IVPMethods/*.cpp")

# Optimize submodule
file(GLOB pybind_LielabOptimizeCPP "lielab/optimize/*.cpp")
file(GLOB pybind_LielabOptimizeExtremizationCPP "lielab/optimize/Extremization/*.cpp")
file(GLOB pybind_LielabOptimizeLeastSquaresCPP "lielab/optimize/LeastSquares/*.cpp")
file(GLOB pybind_LielabOptimizeRootSolvingCPP "lielab/optimize/RootSolving/*.cpp")

# Testing submodule
file(GLOB pybind_LielabTestingCPP "lielab/testing/*.cpp")

# utils submodule
file(GLOB pybind_LielabUtilsCPP "lielab/utils/*.cpp")

pybind11_add_module(${APP_PYD} pybind_Lielab.cpp
    ${pybind_LielabSubmodulesCPP}
    ${pybind_LielabDomainCPP}
    ${pybind_LielabDomainLieAlgebrasCPP}
    ${pybind_LielabDomainLieGroupsCPP}
    ${pybind_LielabDomainSmoothManifoldsCPP}
    ${pybind_LielabFunctionsCPP}
    ${pybind_LielabIntegrateCPP}
    ${pybind_LielabIntegrateIVPMethodsCPP}
    ${pybind_LielabOptimizeCPP}
    ${pybind_LielabOptimizeExtremizationCPP}
    ${pybind_LielabOptimizeLeastSquaresCPP}
    ${pybind_LielabOptimizeRootSolvingCPP}
    ${pybind_LielabTestingCPP}
    ${pybind_LielabUtilsCPP})

if (LIELAB_WITH_ASSERTIONS)
    target_compile_definitions(${APP_PYD} PRIVATE "LIELAB_INCLUDE_ASSERTS")
endif ()

# Link dependencies
target_link_libraries(${APP_PYD} PUBLIC Eigen3::Eigen)
target_link_libraries(${APP_PYD} PUBLIC Lielab::Lielab)

# install(TARGETS ${APP_PYD} DESTINATION .)

# Move the output binary to the Python subfolder
add_custom_command( TARGET ${APP_PYD} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "$<TARGET_FILE:cppLielab>" "${CMAKE_CURRENT_LIST_DIR}/lielab/" )
