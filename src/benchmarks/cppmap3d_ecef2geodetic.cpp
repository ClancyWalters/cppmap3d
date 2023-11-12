#include <doctest/doctest.h>
#include <nanobench.h>
#include "../cppmap3d.hh"

TEST_CASE("ecef2geodetic_benchmark") {
    double lat = 0, lon = 0, alt = 0;

    const auto ELL = cppmap3d::Ellipsoid::WGS84;
    const auto A = cppmap3d::internal::getMajor(ELL);
    const auto B = cppmap3d::internal::getMinor(ELL);
    // clang-format off
    const std::vector<std::vector<std::vector<double>>> data_container {
        {{A, 0, 0}, {0, 0, 0}},
        {{A - 1, 0, 0}, {0, 0, -1}},
        {{A + 1, 0, 0}, {0, 0, 1}},
        {{0.1 * A, 0, 0}, {0, 0, -0.9 * A}},
        //{{0.001 * A, 0, 0}, {0, 0, -0.999 * A}}, test case fails sanity check for height above earth for olson
        {{0, A, 0}, {0, 90, 0}},
        {{0, A - 1, 0}, {0, 90, -1}},
        {{0, A + 1, 0}, {0, 90, 1}},
        {{0, 0.1 * A, 0}, {0, 90, -0.9 * A}},
        //{{0, 0.001 * A, 0}, {0, 90, -0.999 * A}}, test case fails sanity check for height above earth for olson
        {{0, 0, B}, {90, 0, 0}},
        {{0, 0, B + 1}, {90, 0, 1}},
        {{0, 0, B - 1}, {90, 0, -1}},
        {{0, 0, 0.1 * B}, {90, 0, -0.9 * B}},
        //{{0, 0, 0.001 * B}, {90, 0, -0.999 * B}}, test case fails sanity check for height above earth for olson
        {{0, 0, B - 1}, {89.999999, 0, -1}},
        {{0, 0, B - 1}, {89.99999, 0, -1}},
        {{0, 0, -B + 1}, {-90, 0, -1}},
        {{0, 0, -B + 1}, {-89.999999, 0, -1}},
        {{0, 0, -B + 1}, {-89.99999, 0, -1}},
        {{-A + 1, 0, 0}, {0, 180, -1}},
    };
    // clang-format on

    ankerl::nanobench::Bench().minEpochIterations(10000).run(
        "ecef2geodetic_olson",
        [&]() {
            for (const auto& data : data_container) {
                cppmap3d::internal::ecef2geodetic_olson(
                    data[0][0],
                    data[0][1],
                    data[0][2],
                    lat,
                    lon,
                    alt
                );
            }
        }
    );

    ankerl::nanobench::Bench().minEpochIterations(10000).run(
        "ecef2geodetic_you",
        [&]() {
            for (const auto& data : data_container) {
                cppmap3d::ecef2geodetic(
                    data[0][0],
                    data[0][1],
                    data[0][2],
                    lat,
                    lon,
                    alt
                );
            }
        }
    );
}