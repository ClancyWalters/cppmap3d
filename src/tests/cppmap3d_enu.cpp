#include "doctest/doctest.h"

#include <stdexcept>
#include <vector>
#include "../cppmap3d.hh"
#include "cppmap3d_test_util.hh"

TEST_CASE("ecef2enu") {
    const auto ELL = cppmap3d::Ellipsoid::WGS84;
    const auto A = cppmap3d::internal::getMajor(ELL);

    // clang-format off
    std::vector<std::vector<std::vector<double>>> data_container {
        {{0, 0, 0}, {0, 0, 0}, {A, 0, 0}}, 
        {{0, 0, 1000}, {0, 0, 0}, {A + 1000, 0, 0}},
    };
    // clang-format on

    for (const auto& data : data_container) {
        CAPTURE(data);

        double x, y, z;
        cppmap3d::enu2ecef(
            data[0][0],
            data[0][1],
            data[0][2],
            radians(data[1][0]),
            radians(data[1][1]),
            data[1][2],
            x,
            y,
            z
        );

        CHECK(x == doctest::Approx(data[2][0]));
        CHECK(y == doctest::Approx(data[2][1]));
        CHECK(z == doctest::Approx(data[2][2]));

        double e, n, u;
        cppmap3d::ecef2enu(
            x,
            y,
            z,
            radians(data[1][0]),
            radians(data[1][1]),
            data[1][2],
            e,
            n,
            u
        );

        CHECK(e == doctest::Approx(data[0][0]));
        CHECK(n == doctest::Approx(data[0][1]));
        CHECK(u == doctest::Approx(data[0][2]));
    }
}