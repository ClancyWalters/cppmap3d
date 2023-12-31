find_package( doctest REQUIRED )
enable_testing()

add_executable( cppmap3d_tests "main.cpp"
    "../cppmap3d.hh"

    "cppmap3d_tests.cpp"
    "cppmap3d_aer.cpp"
    "cppmap3d_enu.cpp"
    "cppmap3d_geodetic.cpp"
    "cppmap3d_ellipsoid.cpp"
    "cppmap3d_ned.cpp"
)

set_property( TARGET cppmap3d_tests PROPERTY CXX_STANDARD 11 )

target_link_libraries(
    cppmap3d_tests
    PRIVATE
    doctest::doctest
)

include( doctest )
doctest_discover_tests( cppmap3d_tests WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )