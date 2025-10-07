cmake_minimum_required(VERSION 3.23)
set(CMAKE_CXX_STANDARD 20)

project("Lielab"
        VERSION 0.5.0
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

################################
# GLobal project options
################################

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # TODO:
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic)
    add_compile_options(-Wno-unused-parameter) # Ignore unused parameters
    add_compile_options(-fPIC)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # TODO:
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(/W4)
    add_compile_options(/wd4100) # Ignore unused parameters
    add_compile_options(/wd4037) # Ignore std::complex<int>
    add_compile_options(/w44388) # Include signed/unsigned comparison warning
endif()

option(LIELAB_INSTALL_LIBRARY "Install Lielab into default locations." TRUE)

option(LIELAB_BUILD_PYTHON "Build Lielab python module" FALSE)
# set(PYTHON_EXECUTABLE "path to python executable") # Define this to explicitly tell pybind11 which Python to use
option(LIELAB_WITH_COVERAGE "Include test coverage information" FALSE)
option(LIELAB_BUILD_TESTS "Build Lielab tests" FALSE)
# option(LIELAB_BUILD_EXAMPLES "Build Lielab examples" FALSE)



################################
# Make the main project
################################
include(GNUInstallDirs)

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
# file(GLOB LielabIntegrateBVPMethodsCPP "Lielab/integrate/BVPMethods/*.cpp")

# Utils submodule
file(GLOB LielabUtilsCPP "Lielab/utils/*.cpp")

add_library(Lielab STATIC ${LielabDomainCPP}
                          ${LielabDomainLieAlgebrasCPP}
                          ${LielabDomainLieGroupsCPP}
                          ${LielabDomainSmoothManifoldsCPP}
                          ${LielabFunctionsCPP}
                          ${LielabIntegrateCPP}
                          ${LielabIntegrateIVPMethodsCPP}
                        #   ${LielabIntegrateBVPMethodsCPP}
                          ${LielabUtilsCPP})

target_link_libraries(Lielab PRIVATE Eigen3::Eigen)

add_library(Lielab::Lielab ALIAS Lielab)

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


###
# add sources so target shows up in IDE
# file(GLOB Lielab_main_files "Lielab.hpp" "Lielab/*.hpp")
# file(GLOB Lielab_constants_files "Lielab/constants/.hpp")
# file(GLOB Lielab_domain_liealgebras_files "Lielab/domain/liealgebras/*.hpp")
# file(GLOB Lielab_domain_liegroups_files "Lielab/domain/liegroups/*.hpp")
# file(GLOB Lielab_domain_operators_files "Lielab/domain/operators/*.hpp")
# file(GLOB Lielab_domain_files "Lielab/domain/*.hpp")
# file(GLOB Lielab_functions_files "Lielab/functions/*.hpp")
# file(GLOB Lielab_topos_functions_files "Lielab/topos/functions/*.hpp")
# file(GLOB Lielab_topos_files "Lielab/topos/*.hpp")
# file(GLOB Lielab_utils_files "Lielab/utils/*.hpp")
# target_sources(Lielab PUBLIC ${Lielab_main_files})
# target_sources(Lielab PUBLIC ${Lielab_constants_files})
# target_sources(Lielab PUBLIC ${Lielab_domain_liealgebras_files})
# target_sources(Lielab PUBLIC ${Lielab_domain_liegroups_files})
# target_sources(Lielab PUBLIC ${Lielab_domain_operators_files})
# target_sources(Lielab PUBLIC ${Lielab_domain_files})
# target_sources(Lielab PUBLIC ${Lielab_functions_files})
# target_sources(Lielab PUBLIC ${Lielab_topos_functions_files})
# target_sources(Lielab PUBLIC ${Lielab_topos_files})
# target_sources(Lielab PUBLIC ${Lielab_utils_files})

# source_group(TREE "${CMAKE_CURRENT_LIST_DIR}/Lielab" FILES ${Lielab_main_files})
# source_group(TREE "${CMAKE_CURRENT_LIST_DIR}/lielab/src" FILES 
#                                ${Lielab_domain_files} 
#                                ${Lielab_domain_liealgebras_files}
#                                ${Lielab_domain_liegroups_files}
#                                ${Lielab_domain_ops_files}
#                                ${Lielab_functions_files}
#                                ${Lielab_topos_files}
#                                ${Lielab_topos_functions_files})

# Set C++20
target_compile_features(Lielab INTERFACE cxx_std_20)

# Let other CMake targets find Lielab with Lielab::Lielab
export(TARGETS Lielab NAMESPACE Lielab:: FILE LielabTargets.cmake)


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

################################
# Make the examples
################################

# if (LIELAB_BUILD_EXAMPLES)
#     add_subdirectory(examples)
# endif()
