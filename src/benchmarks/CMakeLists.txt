find_package( nanobench REQUIRED )

add_executable( cppmap3d_benchmarks
    "main.cpp"
    "../cppmap3d.hh"
    "nanobench.cpp"
    "cppmap3d_ecef2geodetic.cpp"
)

target_link_libraries(
    cppmap3d_benchmarks
    PRIVATE
    nanobench::nanobench
)