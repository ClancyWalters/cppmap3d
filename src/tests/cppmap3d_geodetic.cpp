#include "doctest/doctest.h"

#include <cmath>
#include <stdexcept>
#include <vector>
#include "../cppmap3d.hh"
#include "cppmap3d_test_util.hh"

static const auto ELL = cppmap3d::Ellipsoid::WGS84;
static const auto A = cppmap3d::internal::getMajor(ELL);
static const auto F = cppmap3d::internal::getFlattening(ELL);
static const auto B = cppmap3d::internal::getMinor(ELL);
static const std::vector<std::vector<std::vector<double>>> xyzlla{
    {{ A, 0, 0 },          { 0, 0, 0 }          },
    { { A - 1, 0, 0 },     { 0, 0, -1 }         },
    { { A + 1, 0, 0 },     { 0, 0, 1 }          },
    { { 0.1 * A, 0, 0 },   { 0, 0, -0.9 * A }   },
    { { 0.001 * A, 0, 0 }, { 0, 0, -0.999 * A } },
    { { 0, A, 0 },         { 0, 90, 0 }         },
    { { 0, A - 1, 0 },     { 0, 90, -1 }        },
    { { 0, A + 1, 0 },     { 0, 90, 1 }         },
    { { 0, 0.1 * A, 0 },   { 0, 90, -0.9 * A }  },
    { { 0, 0.001 * A, 0 }, { 0, 90, -0.999 * A }},
    { { 0, 0, B },         { 90, 0, 0 }         },
    { { 0, 0, B + 1 },     { 90, 0, 1 }         },
    { { 0, 0, B - 1 },     { 90, 0, -1 }        },
    { { 0, 0, 0.1 * B },   { 90, 0, -0.9 * B }  },
    { { 0, 0, 0.001 * B }, { 90, 0, -0.999 * B }},
    { { 0, 0, B - 1 },     { 89.999999, 0, -1 } },
    { { 0, 0, B - 1 },     { 89.99999, 0, -1 }  },
    { { 0, 0, -B + 1 },    { -90, 0, -1 }       },
    { { 0, 0, -B + 1 },    { -89.999999, 0, -1 }},
    { { 0, 0, -B + 1 },    { -89.99999, 0, -1 } },
    { { -A + 1, 0, 0 },    { 0, 180, -1 }       },
};
static const std::vector<std::vector<std::vector<double>>> llaxyz{
    {{ 0, 0, -1 },    { A - 1, 0, 0 } },
    { { 0, 90, -1 },  { 0, A - 1, 0 } },
    { { 0, -90, -1 }, { 0, -A + 1, 0 }},
    { { 90, 0, -1 },  { 0, 0, B - 1 } },
    { { 90, 15, -1 }, { 0, 0, B - 1 } },
    { { -90, 0, -1 }, { 0, 0, -B + 1 }},
};
static const std::vector<double> lla = { radians(42), radians(-82), 200 };
static const std::vector<double> xyz = { 660675.2518247,
                                         -4700948.68316,
                                         4245737.66222 };

TEST_CASE("ecef2geodetic") {
    for (const auto& data : xyzlla) {
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

TEST_CASE("geodetic2ecef") {
    for (const auto& data : llaxyz) {
        CAPTURE(data);

        double x, y, z;
        cppmap3d::geodetic2ecef(
            radians(data[0][0]),
            radians(data[0][1]),
            data[0][2],
            x,
            y,
            z
        );

        CHECK(x == doctest::Approx(data[1][0]));
        CHECK(y == doctest::Approx(data[1][1]));
        CHECK(z == doctest::Approx(data[1][2]));
    }
}

TEST_CASE("ecef2geodetic_inside") {
    double x, y, z;
    double lat, lon, alt;

    CAPTURE(lla);

    cppmap3d::geodetic2ecef(lla[0], lla[1], lla[2], x, y, z);
    cppmap3d::ecef2geodetic(x, y, z, lat, lon, alt);
    CHECK(lat == doctest::Approx(lla[0]));
    CHECK(lon == doctest::Approx(lla[1]));
    CHECK(alt == doctest::Approx(lla[2]));

    cppmap3d::geodetic2ecef(lla[0], lla[1], lla[2], x, y, z);
    cppmap3d::internal::ecef2geodetic_olson(x, y, z, lat, lon, alt);
    CHECK(lat == doctest::Approx(lla[0]));
    CHECK(lon == doctest::Approx(lla[1]));
    CHECK(alt == doctest::Approx(lla[2]));

    cppmap3d::geodetic2ecef(lla[0], lla[1], -lla[2], x, y, z);
    cppmap3d::ecef2geodetic(x, y, z, lat, lon, alt);
    CHECK(lat == doctest::Approx(lla[0]));
    CHECK(lon == doctest::Approx(lla[1]));
    CHECK(alt == doctest::Approx(-lla[2]));

    cppmap3d::geodetic2ecef(lla[0], lla[1], -lla[2], x, y, z);
    cppmap3d::internal::ecef2geodetic_olson(x, y, z, lat, lon, alt);
    CHECK(lat == doctest::Approx(lla[0]));
    CHECK(lon == doctest::Approx(lla[1]));
    CHECK(alt == doctest::Approx(-lla[2]));
}

TEST_CASE("ecef") {
    double x, y, z;
    cppmap3d::geodetic2ecef(lla[0], lla[1], lla[2], x, y, z);

    CHECK(x == doctest::Approx(xyz[0]));
    CHECK(y == doctest::Approx(xyz[1]));
    CHECK(z == doctest::Approx(xyz[2]));

    CHECK_THROWS_AS(
        cppmap3d::geodetic2ecef(-100, lla[1], lla[2], x, y, z),
        std::domain_error
    );
}

TEST_CASE("aer2geodetic") {
    // clang-format off
    std::vector<std::vector<std::vector<double>>> data_container {
        {{33, 77, 1000}, {42.0016981935, -81.99852, 1174.374035}, {42, -82, 200}},
        {{0, 90, 10000}, {0, 0, 10000}, {0, 0, 0}},
    };
    // clang-format on

    for (const auto& data : data_container) {
        CAPTURE(data);

        double lat, lon, alt;
        cppmap3d::aer2geodetic(
            radians(data[0][0]),
            radians(data[0][1]),
            data[0][2],
            radians(data[2][0]),
            radians(data[2][1]),
            data[2][2],
            lat,
            lon,
            alt
        );

        CHECK(lat == doctest::Approx(radians(data[1][0])));
        CHECK(lon == doctest::Approx(radians(data[1][1])));
        CHECK(alt == doctest::Approx(data[1][2]));

        CHECK_THROWS_AS(
            cppmap3d::aer2geodetic(
                radians(data[0][0]),
                radians(data[0][1]),
                -1,
                radians(data[2][0]),
                radians(data[2][1]),
                data[2][2],
                lat,
                lon,
                alt
            ),
            std::domain_error
        );

        double az, el, range;
        cppmap3d::geodetic2aer(
            radians(data[1][0]),
            radians(data[1][1]),
            data[1][2],
            radians(data[2][0]),
            radians(data[2][1]),
            data[2][2],
            az,
            el,
            range
        );
        CHECK(az == doctest::Approx(radians(data[0][0])).epsilon(1e-3));
        CHECK(el == doctest::Approx(radians(data[0][1])).epsilon(1e-3));
        CHECK(range == doctest::Approx(data[0][2]).epsilon(1e-3));
    }
}
