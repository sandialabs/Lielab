cmake_minimum_required(VERSION 3.23)
set(CMAKE_CXX_STANDARD 20)

project("Lielab"
        VERSION 0.5.2
        LANGUAGES CXX)

################################
# Import dependencies
################################

# Manually include dependencies
# add_subdirectory(include/Catch2)
# add_subdirectory(include/Eigen3)
# add_subdirectory(include/pybind11)

if (NOT DEFINED Eigen3_FOUND)
    find_package(Eigen3 REQUIRED)
endif()

if (NOT DEFINED fmt_FOUND)
    find_package(fmt REQUIRED)
endif()

################################
# GLobal project options
################################

# Root "project" options for when the end product is Lielab as a repo.
option(LIELAB_BUILD_PYTHON "Build Lielab python module" FALSE)
# set(PYTHON_EXECUTABLE "path to python executable") # Define this to explicitly tell pybind11 which Python to use
option(LIELAB_WITH_COVERAGE "Include test coverage information" FALSE)
option(LIELAB_BUILD_TESTS "Build Lielab tests" FALSE)
option(LIELAB_BUILD_EXAMPLES "Build Lielab examples" FALSE)

# Lielab as a library options
option(LIELAB_INSTALL_LIBRARY "Install Lielab into default locations." TRUE)

# set(LIELAB_SIMD "auto" CACHE STRING "Vectorization option.")
# set_property(CACHE LIELAB_SIMD PROPERTY STRINGS "auto" "sse2" "avx")

# if(LIELAB_SIMD STREQUAL "sse2")
#     set(_LIELAB_SIMD "sse2")
# elseif(LIELAB_SIMD STREQUAL "avx")
#     set(_LIELAB_SIMD "avx")
# # elseif(LIELAB_SIMD STREQUAL "avx2") # TODO: AVX2 and above blows up compilation
# #     set(_LIELAB_SIMD "avx2")
# # elseif(LIELAB_SIMD STREQUAL "avx512")
# #     set(_LIELAB_SIMD "avx512")
# else()
#     set(_LIELAB_SIMD "auto")
# endif()

################################
# Make the main project
################################
include(GNUInstallDirs)

# Options meta library
add_library(options INTERFACE)
add_library(Lielab::options ALIAS options)

if (LIELAB_WITH_ASSERTIONS)
    target_compile_definitions(options INTERFACE "LIELAB_INCLUDE_ASSERTS")
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    # Warnings
    target_compile_options(options INTERFACE "-Wall")
    target_compile_options(options INTERFACE "-Wextra")
    target_compile_options(options INTERFACE "-Wpedantic")

    target_compile_options(options INTERFACE "-Wno-unused-parameter") # Ignore unused parameters
    target_compile_options(options INTERFACE "-Wno-deprecated-declarations") # Ignore deprecation warning for whenever Eigen::placeholders::last and all are used.

    # Optimization
    target_compile_options(options INTERFACE "-O3")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # TODO:

    # target_compile_options(Lielab PRIVATE -march=native)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Warnings
    target_compile_options(options INTERFACE "-Wall")
    target_compile_options(options INTERFACE "-Wextra")
    target_compile_options(options INTERFACE "-Wpedantic")

    target_compile_options(options INTERFACE "-Wno-unused-parameter") # Ignore unused parameters
    target_compile_options(options INTERFACE "-Wno-deprecated-declarations") # Ignore deprecation warning for whenever Eigen::placeholders::last and all are used.

    # Optimization
    target_compile_options(options INTERFACE "-O3")

    # target_compile_options(options INTERFACE -march=native)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # TODO:
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Warnings
    target_compile_options(options INTERFACE "/W4")
    target_compile_options(options INTERFACE "/wd4100") # Ignore unused parameters.
    target_compile_options(options INTERFACE "/wd4996") # Ignore deprecation warning for whenever Eigen::placeholders::last and ::all are used.
    target_compile_options(options INTERFACE "/wd4037") # Ignore std::complex<int>

    target_compile_options(options INTERFACE "/w44388") # Include signed/unsigned comparison warning

    # Things needed to prevent compilation failures
    target_compile_options(options INTERFACE "/bigobj") # Compilation sometimes fails due to object size

    # Optimization
    target_compile_options(options INTERFACE "$<$<CONFIG:RELEASE>:/O2>") # Optimization. Should be default already.
    target_compile_options(options INTERFACE "$<$<CONFIG:RELEASE>:/Oi>") # Inline intrinsics.
    target_compile_options(options INTERFACE "$<$<CONFIG:RELEASE>:/Qpar>") # Parallelize loops.
    target_compile_options(options INTERFACE "$<$<CONFIG:RELEASE>:/fp:precise>") # Precise floating point model.
    target_compile_options(options INTERFACE "$<$<CONFIG:RELEASE>:/fp:except->") # Ignore precise exceptions.

    # Vectorization
    # Leave this off. It doesnt help performance much but it bloats compile time.
    # It is helpful though since it can ferret out some other bugs.
    # target_compile_options(options INTERFACE "/arch:AVX")
    # target_compile_options(Eigen3::Eigen INTERFACE "/arch:AVX")
endif()

# The rest of Lielab

# Root
file(GLOB LielabRootCPP "Lielab/*.cpp")

# Domain submodule
file(GLOB LielabDomainCPP "Lielab/domain/*.cpp")
file(GLOB LielabDomainLieAlgebrasCPP "Lielab/domain/liealgebras/*.cpp")
file(GLOB LielabDomainLieGroupsCPP "Lielab/domain/liegroups/*.cpp")
file(GLOB LielabDomainSmoothManifoldsCPP "Lielab/domain/smoothmanifolds/*.cpp")

# Functions submodule
file(GLOB LielabFunctionsCPP "Lielab/functions/*.cpp")

# Integrate submodule
file(GLOB LielabIntegrateCPP "Lielab/integrate/*.cpp")
file(GLOB LielabIntegrateIVPMethodsCPP "Lielab/integrate/IVPMethods/*.cpp")

# Optimize submodule
file(GLOB LielabOptimizeCPP "Lielab/optimize/*.cpp")
file(GLOB LielabOptimizeExtremizationCPP "Lielab/optimize/Extremization/*.cpp")
file(GLOB LielabOptimizeRootSolvingCPP "Lielab/optimize/RootSolving/*.cpp")

# Testing submodule
file(GLOB LielabTestingCPP "Lielab/testing/*.cpp")

# Utils submodule
file(GLOB LielabUtilsCPP "Lielab/utils/*.cpp")

add_library(Lielab STATIC ${LielabRootCPP}
                          ${LielabDomainCPP}
                          ${LielabDomainLieAlgebrasCPP}
                          ${LielabDomainLieGroupsCPP}
                          ${LielabDomainSmoothManifoldsCPP}
                          ${LielabFunctionsCPP}
                          ${LielabIntegrateCPP}
                          ${LielabIntegrateIVPMethodsCPP}
                          ${LielabOptimizeCPP}
                          ${LielabOptimizeExtremizationCPP}
                          ${LielabOptimizeRootSolvingCPP}
                          ${LielabTestingCPP}
                          ${LielabUtilsCPP})
add_library(Lielab::Lielab ALIAS Lielab)

set_target_properties(Lielab PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

target_link_libraries(Lielab PUBLIC Eigen3::Eigen)
target_link_libraries(Lielab PRIVATE fmt::fmt)
target_link_libraries(Lielab PRIVATE Lielab::options)

# Link options
# if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#     # TODO:
# elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#     # TODO:
# elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
#     # TODO:
# elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
#     # TODO:
# endif()

target_include_directories(Lielab PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Lielab>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if (LIELAB_INSTALL_LIBRARY)
    # TODO: Are these needed?
    # PUBLIC_HEADER DESTINATION include COMPONENT Development
    # BUNDLE DESTINATION bin COMPONENT Runtime
    install(TARGETS Lielab
            EXPORT LielabTargets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    include(CMakePackageConfigHelpers)
    # configure_package_config_file(
    #     "${PROJECT_SOURCE_DIR}/cmake/LielabConfig.cmake.in"
    #     "${PROJECT_BINARY_DIR}/LielabConfig.cmake"
    #     INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Lielab/cmake
    # )

    # install(EXPORT LielabTargets
    #         FILE LielabTargets.cmake
    #         NAMESPACE Lielab::
    #         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Lielab/cmake)

    # install(FILES "${PROJECT_BINARY_DIR}/LielabConfigVersion.cmake"
    #               "${PROJECT_BINARY_DIR}/LielabConfig.cmake"
    #         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Lielab/cmake)

    # Copy the extra files in the root
    install(FILES Lielab.hpp
            DESTINATION include)
    install(FILES LICENSE
            DESTINATION licenses)

    # Copy the whole "lielab" directory
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/Lielab
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
                       PATTERN "*.tpp")
endif()

# Set C++20
target_compile_features(Lielab INTERFACE cxx_std_20)

target_precompile_headers(Lielab PRIVATE <Lielab/pch.hpp>)

# Let other CMake targets find Lielab with Lielab::Lielab
export(TARGETS options Lielab NAMESPACE Lielab:: FILE LielabTargets.cmake)
set(Lielab_FOUND TRUE)


################################
# Make the tests project
################################
if (LIELAB_BUILD_TESTS)
    if (LIELAB_WITH_COVERAGE AND (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"))
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage -fcondition-coverage")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 --coverage -fcondition-coverage")
        target_link_libraries(Lielab PRIVATE gcov)
    endif ()
    add_subdirectory(tests)
endif()

################################
# Make the Python library
################################
if (LIELAB_BUILD_PYTHON)
    add_subdirectory(python)
endif()
