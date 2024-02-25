project(pyLielab)

# Find python
set(CMAKE_CXX_STANDARD 20)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Add the pybind11 cpp file
SET( APP_EXE cppLielab )
pybind11_add_module(${APP_EXE} cppLielab.cpp)

if(MSVC)
    target_compile_definitions(cppLielab PUBLIC _USE_MATH_DEFINES)
    target_compile_options(cppLielab PUBLIC -bigobj)
endif(MSVC)

target_link_libraries(${APP_EXE} PUBLIC Eigen3::Eigen)
target_link_libraries(${APP_EXE} PUBLIC Lielab::Lielab)
    
add_custom_command( TARGET cppLielab POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "$<TARGET_FILE:cppLielab>" "${CMAKE_CURRENT_LIST_DIR}/lielab/" )
