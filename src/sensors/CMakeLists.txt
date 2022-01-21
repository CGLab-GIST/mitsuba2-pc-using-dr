set(MTS_PLUGIN_PREFIX "sensors")

add_plugin(perspective     perspective.cpp)
add_plugin(perspective_ex     perspective_ex.cpp)


# Register the test directory
add_tests(${CMAKE_CURRENT_SOURCE_DIR}/tests)
