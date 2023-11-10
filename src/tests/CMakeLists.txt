add_executable( cppmap3d_tests "main.cpp" "../cppmap3d.hh" "cppmap3d_tests.cpp" )

target_link_libraries(
    cppmap3d_tests
    PRIVATE
    doctest::doctest
)

include( doctest )
doctest_discover_tests( cppmap3d_tests WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )