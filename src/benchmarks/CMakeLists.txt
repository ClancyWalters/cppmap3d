include(FetchContent)

FetchContent_Declare(
    nanobench
    GIT_REPOSITORY https://github.com/martinus/nanobench
    GIT_TAG v4.3.11
)
FetchContent_MakeAvailable(nanobench)

FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest
    GIT_TAG v2.4.11
)
FetchContent_MakeAvailable(doctest)

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
    doctest::doctest
)