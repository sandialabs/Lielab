project(tests LANGUAGES CXX)
# include(CTest)

enable_testing()

if (NOT DEFINED Catch2_FOUND)
    find_package(Catch2 REQUIRED)
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # TODO:
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # TODO:
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # TODO:
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(/wd4245) # Ignore lazy casting int to size_t
    add_compile_options(/wd4037) # Ignore std::complex<int>
    add_compile_options(/wd4996) # Ignore lazy casting size_t to int
endif()

# Add executable and cpp files
file(GLOB TestLielabDomainLieAlgebrasCPP "domain/liealgebras/*.cpp")
file(GLOB TestLielabDomainLieGroupsCPP "domain/liegroups/*.cpp")
file(GLOB TestLielabDomainSmoothManifoldsCPP "domain/smoothmanifolds/*.cpp")
file(GLOB TestLielabFunctionsCPP "functions/*.cpp")
file(GLOB TestLielabIntegrateCPP "integrate/*.cpp")
file(GLOB TestLielabOptimizeCPP "optimize/*.cpp")
file(GLOB TestLielabTestingCPP "testing/*.cpp")
file(GLOB TestLielabUtilsCPP "utils/*.cpp")

# SET(test_Lielab test_Lielab)
ADD_EXECUTABLE(test_Lielab test_Lielab.cpp
        ${TestLielabDomainLieAlgebrasCPP}
        ${TestLielabDomainLieGroupsCPP}
        ${TestLielabDomainSmoothManifoldsCPP}
        ${TestLielabFunctionsCPP}
        ${TestLielabIntegrateCPP}
        ${TestLielabOptimizeCPP}
        ${TestLielabTestingCPP}
        ${TestLielabUtilsCPP})

# Link dependencies
target_link_libraries(test_Lielab PRIVATE Catch2::Catch2)
target_link_libraries(test_Lielab PRIVATE Lielab::Lielab)
target_link_libraries(test_Lielab PRIVATE Eigen3::Eigen)

# Link options
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # TODO:
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # TODO:
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # TODO:
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set_target_properties(test_Lielab PROPERTIES LINK_FLAGS "/PROFILE")
endif()
