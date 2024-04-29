cmake_minimum_required(VERSION 3.23)

project("Lielab"
        VERSION 0.3.0
        LANGUAGES CXX)

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Is top level, compile everything.
    # Is "PROJECT_IS_TOP_LEVEL" the same thing?
    set(IS_TOPLEVEL_PROJECT TRUE)
else()
    # Not top level. Imported into another project.
    set(IS_TOPLEVEL_PROJECT FALSE)
endif()

option(LIELAB_INSTALL_LIBRARY "Install Lielab into default locations." ${IS_TOPLEVEL_PROJECT})
option(LIELAB_BUILD_TESTS "Build Lielab tests" ${IS_TOPLEVEL_PROJECT})
option(LIELAB_BUILD_PYTHON "Build Lielab python module" ${IS_TOPLEVEL_PROJECT})

################################
# Make the main project
#
# Always build this regardless of how this CMakeLists is called.
################################
include(GNUInstallDirs)

add_library(Lielab INTERFACE)
add_library(Lielab::Lielab ALIAS Lielab)

target_include_directories(Lielab INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Lielab>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)


# Uncomment these if you want to manually include dependencies and build with CMake.
# add_subdirectory(include/Catch2)
# add_subdirectory(include/Eigen3)
# add_subdirectory(include/pybind11)

# Eigen3 required for the header only pure C++ version
# if (NOT DEFINED Eigen3_FOUND)
#     find_package(Eigen3 REQUIRED)
# endif()

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
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/Lielab DESTINATION include)
endif()


###
# add sources so target shows up in IDE
# file(GLOB Lielab_main_files "lielab/*.hpp" "lielab/lielab" "lielab/topos" "lielab/functions" "lielab/domain")
# file(GLOB Lielab_domain_files "lielab/src/domain/*.hpp")
# file(GLOB Lielab_domain_liealgebras_files "lielab/src/domain/liealgebras/*.hpp")
# file(GLOB Lielab_domain_liegroups_files "lielab/src/domain/liegroups/*.hpp")
# file(GLOB Lielab_domain_ops_files "lielab/src/domain/ops/*.hpp")
# file(GLOB Lielab_functions_files "lielab/src/functions/*.hpp")
# file(GLOB Lielab_topos_files "lielab/src/topos/*.hpp")
# file(GLOB Lielab_topos_functions_files "lielab/src/topos/functions/*.hpp")
# target_sources(Lielab PRIVATE ${Lielab_main_files})
# target_sources(Lielab PRIVATE ${Lielab_domain_files})
# target_sources(Lielab PRIVATE ${Lielab_domain_liealgebras_files})
# target_sources(Lielab PRIVATE ${Lielab_domain_liegroups_files})
# target_sources(Lielab PRIVATE ${Lielab_domain_ops_files})
# target_sources(Lielab PRIVATE ${Lielab_functions_files})
# target_sources(Lielab PRIVATE ${Lielab_topos_files})
# target_sources(Lielab PRIVATE ${Lielab_topos_functions_files})

# source_group(TREE "${CMAKE_CURRENT_LIST_DIR}/lielab" FILES ${Lielab_main_files})
# source_group(TREE "${CMAKE_CURRENT_LIST_DIR}/lielab/src" FILES 
#                                ${Lielab_domain_files} 
#                                ${Lielab_domain_liealgebras_files}
#                                ${Lielab_domain_liegroups_files}
#                                ${Lielab_domain_ops_files}
#                                ${Lielab_functions_files}
#                                ${Lielab_topos_files}
#                                ${Lielab_topos_functions_files})

# if(MSVC)
#     target_compile_definitions(Lielab INTERFACE "_USE_MATH_DEFINES")
#     target_compile_options(Lielab INTERFACE "-bigobj")
#     target_compile_options(Lielab INTERFACE "/fp:precise")
#     target_compile_definitions(Lielab INTERFACE "_CRT_SECURE_NO_WARNINGS")
#     target_compile_options(Lielab INTERFACE "/Zc:preprocessor") # use new standard-conforming preprocessor by default
#     # CMake does not set include path correctly for interface targets
#     # get_target_property(EIGEN3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
#     get_target_property(LIELAB_INCLUDE_DIRS Lielab INTERFACE_INCLUDE_DIRECTORIES)
#     foreach(dir ${EIGEN3_INCLUDE_DIRS} ${LIELAB_INCLUDE_DIRS})
#         string(APPEND LIELAB_INCLUDE_DIR_STRING "${dir};")
#     endforeach()
#     set_target_properties(Lielab PROPERTIES VS_GLOBAL_IncludePath "${LIELAB_INCLUDE_DIR_STRING}")
# endif(MSVC)
###

# Set C++20
target_compile_features(Lielab INTERFACE cxx_std_20)

# Let other CMake targets find Lielab with Lielab::Lielab
export(TARGETS Lielab NAMESPACE Lielab:: FILE LielabTargets.cmake)


################################
# Make the tests project
################################
if (LIELAB_BUILD_TESTS)
    if (NOT DEFINED Eigen3_FOUND)
        find_package(Eigen3 REQUIRED)
    endif()
    if (NOT DEFINED Catch2_FOUND)
        find_package(Catch2 REQUIRED)
    endif()
    add_subdirectory(tests)
endif()


################################
# Make the Python library
################################
if (LIELAB_BUILD_PYTHON)
    if (NOT DEFINED Eigen3_FOUND)
        find_package(Eigen3 REQUIRED)
    endif()
    if (NOT DEFINED pybind11_FOUND)
        find_package(pybind11 REQUIRED)
    endif()
    add_subdirectory(python)
endif()

################################
# 32-bit cl.exe will run out of memory building Lielab
################################
if (MSVC AND (NOT CMAKE_VS_PLATFORM_TOOLSET_HOST_ARCHITECTURE MATCHES "x64"))
    message(WARNING "Use /Ax64 and/or /Thost=x64 to make sure 64-bit host toolchain is used. CMAKE_VS_PLATFORM_TOOLSET_HOST_ARCHITECTURE is ${CMAKE_VS_PLATFORM_TOOLSET_HOST_ARCHITECTURE}")
endif()
