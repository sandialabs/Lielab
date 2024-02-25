project(tests)
include(CTest)

# Add executable
SET( APP_EXE Lielab_test )
ADD_EXECUTABLE( ${APP_EXE} test_main.cpp )

# Link dependencies
target_link_libraries(${APP_EXE} PRIVATE Catch2::Catch2WithMain)
target_link_libraries(${APP_EXE} PRIVATE Lielab::Lielab)
target_link_libraries(${APP_EXE} PRIVATE Eigen3::Eigen)