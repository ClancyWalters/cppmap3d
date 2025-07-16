#include "doctest/doctest.h"

#include <stdexcept>
#include <vector>
#include "../cppmap3d.hh"
#include "cppmap3d_test_util.hh"

static const std::vector<double> lla = { radians(42), radians(-82), 200 };
static const std::vector<double> aer = { 33, 70, 1000 };

static const auto ELL = cppmap3d::Ellipsoid::WGS84;
static const auto A = cppmap3d::internal::getMajor(ELL);
static const auto B = cppmap3d::internal::getMinor(ELL);

TEST_CASE("ecef_ned") {
    double enu_e, enu_n, enu_u;
    cppmap3d::aer2enu(aer[0], aer[1], aer[2], enu_e, enu_n, enu_u);

    std::vector<double> ned = { enu_n, enu_e, -enu_u };

    double x, y, z;
    cppmap3d::aer2ecef(aer[0], aer[1], aer[2], lla[0], lla[1], lla[2], x, y, z);

    double n, e, d;
    cppmap3d::ecef2ned(x, y, z, lla[0], lla[1], lla[2], n, e, d);

    CHECK(n == doctest::Approx(ned[0]));
    CHECK(e == doctest::Approx(ned[1]));
    CHECK(d == doctest::Approx(ned[2]));

    double x1, y1, z1;
    cppmap3d::ned2ecef(
        ned[0],
        ned[1],
        ned[2],
        lla[0],
        lla[1],
        lla[2],
        x1,
        y1,
        z1
    );

    CHECK(x == doctest::Approx(x1));
    CHECK(y == doctest::Approx(y1));
    CHECK(z == doctest::Approx(z1));
}

TEST_CASE("ned_geodetic") {
    double lat, lon, alt;
    cppmap3d::aer2geodetic(
        aer[0],
        aer[1],
        aer[2],
        lla[0],
        lla[1],
        lla[2],
        lat,
        lon,
        alt
    );

    double enu_e, enu_n, enu_u;
    cppmap3d::geodetic2enu(
        lat,
        lon,
        alt,
        lla[0],
        lla[1],
        lla[2],
        enu_e,
        enu_n,
        enu_u
    );

    std::vector<double> ned = { enu_n, enu_e, -enu_u };

    double n, e, d;
    cppmap3d::geodetic2ned(lat, lon, alt, lla[0], lla[1], lla[2], n, e, d);

    CHECK(n == doctest::Approx(ned[0]));
    CHECK(e == doctest::Approx(ned[1]));
    CHECK(d == doctest::Approx(ned[2]));

    double lat1, lon1, alt1;
    cppmap3d::enu2geodetic(
        enu_e,
        enu_n,
        enu_u,
        lla[0],
        lla[1],
        lla[2],
        lat1,
        lon1,
        alt1
    );

    CHECK(lat == doctest::Approx(lat1));
    CHECK(lon == doctest::Approx(lon1));
    CHECK(alt == doctest::Approx(alt1));

    cppmap3d::ned2geodetic(
        ned[0],
        ned[1],
        ned[2],
        lla[0],
        lla[1],
        lla[2],
        lat1,
        lon1,
        alt1
    );

    CHECK(lat == doctest::Approx(lat1));
    CHECK(lon == doctest::Approx(lon1));
    CHECK(alt == doctest::Approx(alt1));
}

TEST_CASE("enuv_nedv") {
    double vx = 5; 
    double vy = 3;
    double vz = 2;
    double ve = 5.368859646588048;
    double vn = 3.008520763668120;
    double vu = -0.352347711524077;

    double e, n, u;
    cppmap3d::ecef2enuv(vx, vy, vz, lla[0], lla[1], e, n, u);

    CHECK(ve == doctest::Approx(e));
    CHECK(vn == doctest::Approx(n));
    CHECK(vu == doctest::Approx(u));

    double x, y, z;
    cppmap3d::enu2ecefv(ve, vn, vu, lla[0], lla[1], x, y, z);

    CHECK(vx == doctest::Approx(x));
    CHECK(vy == doctest::Approx(y));
    CHECK(vz == doctest::Approx(z));

    n = 0;
    e = 0;
    double d;
    cppmap3d::ecef2nedv(vx, vy, vz, lla[0], lla[1], n, e, d);

    CHECK(vn == doctest::Approx(n));
    CHECK(ve == doctest::Approx(e));
    CHECK(-vu == doctest::Approx(d));
}