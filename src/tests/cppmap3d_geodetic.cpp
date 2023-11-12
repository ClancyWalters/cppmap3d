#include "doctest/doctest.h"

#include <stdexcept>
#include <vector>
#include "../cppmap3d.hh"
#include "cppmap3d_test_util.hh"

TEST_CASE("ecef2geodetic") {
    const auto ELL = cppmap3d::Ellipsoid::WGS84;
    const auto A = cppmap3d::internal::getMajor(ELL);
    const auto F = cppmap3d::internal::getFlattening(ELL);
    const auto B = cppmap3d::internal::getMinor(ELL);

    // clang-format off
    std::vector<std::vector<std::vector<double>>> data_container {
        {{A, 0, 0}, {0, 0, 0}},
        {{A - 1, 0, 0}, {0, 0, -1}},
        {{A + 1, 0, 0}, {0, 0, 1}},
        {{0.1 * A, 0, 0}, {0, 0, -0.9 * A}},
        {{0.001 * A, 0, 0}, {0, 0, -0.999 * A}},
        {{0, A, 0}, {0, 90, 0}},
        {{0, A - 1, 0}, {0, 90, -1}},
        {{0, A + 1, 0}, {0, 90, 1}},
        {{0, 0.1 * A, 0}, {0, 90, -0.9 * A}},
        {{0, 0.001 * A, 0}, {0, 90, -0.999 * A}},
        {{0, 0, B}, {90, 0, 0}},
        {{0, 0, B + 1}, {90, 0, 1}},
        {{0, 0, B - 1}, {90, 0, -1}},
        {{0, 0, 0.1 * B}, {90, 0, -0.9 * B}},
        {{0, 0, 0.001 * B}, {90, 0, -0.999 * B}},
        {{0, 0, B - 1}, {89.999999, 0, -1}},
        {{0, 0, B - 1}, {89.99999, 0, -1}},
        {{0, 0, -B + 1}, {-90, 0, -1}},
        {{0, 0, -B + 1}, {-89.999999, 0, -1}},
        {{0, 0, -B + 1}, {-89.99999, 0, -1}},
        {{-A + 1, 0, 0}, {0, 180, -1}},
    };
    // clang-format on

    for (const auto& data : data_container) {
        CAPTURE(data);

        double lat, lon, alt;
        cppmap3d::ecef2geodetic(
            data[0][0],
            data[0][1],
            data[0][2],
            lat,
            lon,
            alt
        );

        CHECK(lat == doctest::Approx(radians(data[1][0])));
        CHECK(lon == doctest::Approx(radians(data[1][1])));
        CHECK(alt == doctest::Approx(data[1][2]));
    }
}

TEST_CASE("ecef2geodetic_olson") {
    const auto ELL = cppmap3d::Ellipsoid::WGS84;
    const auto A = cppmap3d::internal::getMajor(ELL);
    const auto F = cppmap3d::internal::getFlattening(ELL);
    const auto B = cppmap3d::internal::getMinor(ELL);

    // clang-format off
    std::vector<std::vector<std::vector<double>>> data_container {
        {{A, 0, 0}, {0, 0, 0}},
        {{A - 1, 0, 0}, {0, 0, -1}},
        {{A + 1, 0, 0}, {0, 0, 1}},
        {{0.1 * A, 0, 0}, {0, 0, -0.9 * A}},
        //{{0.001 * A, 0, 0}, {0, 0, -0.999 * A}}, test case fails sanity check for height above earth
        {{0, A, 0}, {0, 90, 0}},
        {{0, A - 1, 0}, {0, 90, -1}},
        {{0, A + 1, 0}, {0, 90, 1}},
        {{0, 0.1 * A, 0}, {0, 90, -0.9 * A}},
        //{{0, 0.001 * A, 0}, {0, 90, -0.999 * A}}, test case fails sanity check for height above earth
        {{0, 0, B}, {90, 0, 0}},
        {{0, 0, B + 1}, {90, 0, 1}},
        {{0, 0, B - 1}, {90, 0, -1}},
        {{0, 0, 0.1 * B}, {90, 0, -0.9 * B}},
        //{{0, 0, 0.001 * B}, {90, 0, -0.999 * B}}, test case fails sanity check for height above earth
        {{0, 0, B - 1}, {89.999999, 0, -1}},
        {{0, 0, B - 1}, {89.99999, 0, -1}},
        {{0, 0, -B + 1}, {-90, 0, -1}},
        {{0, 0, -B + 1}, {-89.999999, 0, -1}},
        {{0, 0, -B + 1}, {-89.99999, 0, -1}},
        {{-A + 1, 0, 0}, {0, 180, -1}},
    };
    // clang-format on

    for (const auto& data : data_container) {
        CAPTURE(data);

        double lat, lon, alt;
        cppmap3d::internal::ecef2geodetic_olson(
            data[0][0],
            data[0][1],
            data[0][2],
            lat,
            lon,
            alt
        );

        CHECK(lat == doctest::Approx(radians(data[1][0])));
        CHECK(lon == doctest::Approx(radians(data[1][1])));
        CHECK(alt == doctest::Approx(data[1][2]));
    }
}

TEST_CASE("ecef2geodetic_benchmark") {}