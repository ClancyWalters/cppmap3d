#include "doctest/doctest.h"

#include <stdexcept>
#include <vector>
#include "../cppmap3d.hh"
#include "cppmap3d_test_util.hh"

TEST_CASE("aer2ecef") {
    // clang-format off
    std::vector<std::vector<std::vector<double>>> data_container {
        {{ 33, 70, 1000 }, { 42, -82, 200 }, { 660930.2, -4701424.0, 4246579.6 }},
    };
    // clang-format on

    for (const auto& data : data_container) {
        CAPTURE(data);

        double x, y, z;
        cppmap3d::aer2ecef(
            radians(data[0][0]),
            radians(data[0][1]),
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

        CHECK_THROWS_AS(
            cppmap3d::aer2ecef(
                radians(data[0][0]),
                radians(data[0][1]),
                -1,
                radians(data[1][0]),
                radians(data[1][1]),
                data[1][2],
                x,
                y,
                z
            ),
            std::domain_error
        );
    }
}

TEST_CASE("ecef2aer") {
    const auto ELL = cppmap3d::Ellipsoid::WGS84;
    const auto A = cppmap3d::internal::getMajor(ELL);
    const auto B = cppmap3d::internal::getMinor(ELL);
    // clang-format off
    std::vector<std::vector<std::vector<double>>> data_container {
        {{A - 1, 0, 0}, {0, 0, 0}, {0, -90, 1}},
        {{-A + 1, 0, 0}, {0, 180, 0}, {0, -90, 1}},
        {{0, A - 1, 0}, {0, 90, 0}, {0, -90, 1}},
        {{0, -A + 1, 0}, {0, -90, 0}, {0, -90, 1}},
        {{0, 0, B - 1}, {90, 0, 0}, {0, -90, 1}},
        {{0, 0, -B + 1}, {-90, 0, 0}, {0, -90, 1}},
        {{660930.19276, -4701424.22296, 4246579.60463}, {42, -82, 200}, {33, 70, 1000}},
    };
    // clang-format on

    for (const auto& data : data_container) {
        CAPTURE(data);

        double a, e, r;
        cppmap3d::ecef2aer(
            data[0][0],
            data[0][1],
            data[0][2],
            radians(data[1][0]),
            radians(data[1][1]),
            data[1][2],
            a,
            e,
            r
        );

        CHECK(a == doctest::Approx(radians(data[2][0])));
        CHECK(e == doctest::Approx(radians(data[2][1])));
        CHECK(r == doctest::Approx(data[2][2]));
    }
}

TEST_CASE("aer2enu") {
    // clang-format off
    std::vector<std::vector<std::vector<double>>> data_container {
        {{ 33, 70, 1000 }, { 186.2775, 286.8422, 939.6926 }},
    };
    // clang-format on

    for (const auto& data : data_container) {
        CAPTURE(data);

        double e, n, u;
        cppmap3d::aer2enu(
            radians(data[0][0]),
            radians(data[0][1]),
            data[0][2],
            e,
            n,
            u
        );

        CHECK(e == doctest::Approx(data[1][0]));
        CHECK(n == doctest::Approx(data[1][1]));
        CHECK(u == doctest::Approx(data[1][2]));

        CHECK_THROWS_AS(
            cppmap3d::aer2enu(
                radians(data[0][0]),
                radians(data[0][1]),
                -1,
                e,
                n,
                u
            ),
            std::domain_error
        );
    }
}

TEST_CASE("aer2ned") {
    // clang-format off
    std::vector<std::vector<std::vector<double>>> data_container {
        {{ 33, 70, 1000 }, { 286.8422, 186.2775, -939.6926 }},
    };
    // clang-format on

    for (const auto& data : data_container) {
        CAPTURE(data);

        double n, e, d;
        cppmap3d::aer2ned(
            radians(data[0][0]),
            radians(data[0][1]),
            data[0][2],
            n,
            e,
            d
        );

        CHECK(n == doctest::Approx(data[1][0]));
        CHECK(e == doctest::Approx(data[1][1]));
        CHECK(d == doctest::Approx(data[1][2]));

        CHECK_THROWS_AS(
            cppmap3d::aer2ned(
                radians(data[0][0]),
                radians(data[0][1]),
                -1,
                n,
                e,
                d
            ),
            std::domain_error
        );
    }
}