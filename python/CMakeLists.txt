project(pyLielab LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)


# Add the pybind11 cpp file
SET( APP_PYD cppLielab )
pybind11_add_module(${APP_PYD} cppLielab.cpp)

if(MSVC)
    target_compile_definitions(${APP_PYD} PUBLIC _USE_MATH_DEFINES)
    target_compile_options(${APP_PYD} PUBLIC -bigobj)
endif(MSVC)

# Link dependencies
target_link_libraries(${APP_PYD} PUBLIC Eigen3::Eigen)
target_link_libraries(${APP_PYD} PUBLIC Lielab::Lielab)

# Move the output binary to the Python subfolder
add_custom_command( TARGET ${APP_PYD} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "$<TARGET_FILE:cppLielab>" "${CMAKE_CURRENT_LIST_DIR}/lielab/" )
