// Tests ported from https://github.com/gberrante/map_3d

#include "doctest/doctest.h"

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "../cppmap3d.hh"
#include "cppmap3d_test_util.hh"

TEST_CASE("geodetic2ecef") {
    double lat = radians(30.14988205);
    double lon = radians(91.38733072);
    double alt = 4031.0;

    double x, y, z;
    cppmap3d::geodetic2ecef(lat, lon, alt, x, y, z);

    double xref = -1.337281037300386e+05;
    double yref = 5.521796910920261e+06;
    double zref = 3.186776473672415e+06;

    CHECK(std::abs(x - xref) < 1e-3);
    CHECK(std::abs(y - yref) < 1e-3);
    CHECK(std::abs(z - zref) < 1e-3);

    CHECK(x == doctest::Approx(xref));
    CHECK(y == doctest::Approx(yref));
    CHECK(z == doctest::Approx(zref));
}

TEST_CASE("geodetic2aer") {
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;

    double lat = radians(42.002581974253744);
    double lon = radians(-81.997751960067460);
    double alt = 1.139701799575106e+03;

    double azref = radians(32.999999999989740);
    double elref = radians(69.999999999945540);
    double rangeref = 1000.0;

    double a, e, r;
    cppmap3d::geodetic2aer(lat, lon, alt, lat0, lon0, alt0, a, e, r);

    CHECK(std::abs(a - azref) < 1e-3);
    CHECK(std::abs(e - elref) < 1e-3);
    CHECK(std::abs(r - rangeref) < 1e-3);

    CHECK(a == doctest::Approx(azref));
    CHECK(e == doctest::Approx(elref));
    CHECK(r == doctest::Approx(rangeref));
}

TEST_CASE("geodetic2enu") {
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;

    double lat = radians(42.002581974253744);
    double lon = radians(-81.997751960067460);
    double alt = 1.139701799575106e+03;

    double eref = 1.862775208168244e+02;
    double nref = 2.868422278521820e+02;
    double uref = 9.396926207845534e+02;

    double e, n, u;
    cppmap3d::geodetic2enu(lat, lon, alt, lat0, lon0, alt0, e, n, u);

    CHECK(std::abs(e - eref) < 1e-3);
    CHECK(std::abs(n - nref) < 1e-3);
    CHECK(std::abs(u - uref) < 1e-3);

    CHECK(e == doctest::Approx(eref));
    CHECK(n == doctest::Approx(nref));
    CHECK(u == doctest::Approx(uref));
}

TEST_CASE("aer2ecef") {
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;

    double az = radians(33);
    double el = radians(70);
    double range = 1000.0;

    double x, y, z;
    cppmap3d::aer2ecef(az, el, range, lat0, lon0, alt0, x, y, z);

    double xref = 6.609301927610816e+05;
    double yref = -4.701424222957011e+06;
    double zref = 4.246579604632881e+06;

    CHECK(std::abs(x - xref) < 1e-3);
    CHECK(std::abs(y - yref) < 1e-3);
    CHECK(std::abs(z - zref) < 1e-3);

    CHECK(x == doctest::Approx(xref));
    CHECK(y == doctest::Approx(yref));
    CHECK(z == doctest::Approx(zref));
}

TEST_CASE("aer2enu") {
    double az = radians(33);
    double el = radians(70);
    double range = 1000.0;

    double eref = 1.862775208165935e+02;
    double nref = 2.868422278517140e+02;
    double uref = 9.396926207859083e+02;

    double e, n, u;
    cppmap3d::aer2enu(az, el, range, e, n, u);

    CHECK(std::abs(e - eref) < 1e-3);
    CHECK(std::abs(n - nref) < 1e-3);
    CHECK(std::abs(u - uref) < 1e-3);

    CHECK(e == doctest::Approx(eref));
    CHECK(n == doctest::Approx(nref));
    CHECK(u == doctest::Approx(uref));
}

TEST_CASE("aer2geodetic") {
    double lat0 = radians(42);
    double lon0 = radians(-82);
    double alt0 = 200.0;

    double az = radians(32.999999999989740);
    double el = radians(69.999999999945540);
    double range = 1000.0;

    double latref = radians(42.002581974253744);
    double lonref = radians(-81.997751960067460);
    double altref = 1.139701799575106e+03;

    double lat, lon, alt;
    cppmap3d::aer2geodetic(az, el, range, lat0, lon0, alt0, lat, lon, alt);

    CHECK(std::abs(lat - latref) < 1e-8);
    CHECK(std::abs(lon - lonref) < 1e-8);
    CHECK(std::abs(alt - altref) < 1e-8);

    CHECK(lat == doctest::Approx(latref));
    CHECK(lon == doctest::Approx(lonref));
    CHECK(alt == doctest::Approx(altref));
}

TEST_CASE("enu2aer") {
    double e = 1.862775210000000e+02;
    double n = 2.868422200000000e+02;
    double u = 9.396926200000000e+02;

    double azref = radians(33.0);
    double elref = radians(70.0);
    double rangeref = 1000.0;

    double az, el, range;
    cppmap3d::enu2aer(e, n, u, az, el, range);

    CHECK(std::abs(az - azref) < 1e-3);
    CHECK(std::abs(el - elref) < 1e-3);
    CHECK(std::abs(range - rangeref) < 1e-3);

    CHECK(az == doctest::Approx(azref));
    CHECK(el == doctest::Approx(elref));
    CHECK(range == doctest::Approx(rangeref));
}

TEST_CASE("enu2ecef") {
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;

    double e = 1.862775210000000e+02;
    double n = 2.868422200000000e+02;
    double u = 9.396926200000000e+02;

    double xref = 6.609301927610815e+05;
    double yref = -4.701424222957011e+06;
    double zref = 4.246579604632881e+06;

    double x, y, z;
    cppmap3d::enu2ecef(e, n, u, lat0, lon0, alt0, x, y, z);

    CHECK(std::abs(x - xref) < 1e-3);
    CHECK(std::abs(y - yref) < 1e-3);
    CHECK(std::abs(z - zref) < 1e-3);

    CHECK(x == doctest::Approx(xref));
    CHECK(y == doctest::Approx(yref));
    CHECK(z == doctest::Approx(zref));
}

TEST_CASE("enu2geodetic") {
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;

    double e = 0.0;
    double n = 0.0;
    double u = -1.0;

    double latref = radians(41.999999999999993);
    double lonref = radians(-82.0);
    double altref = 1.990000000007368e+02;

    double lat, lon, alt;
    cppmap3d::enu2geodetic(e, n, u, lat0, lon0, alt0, lat, lon, alt);

    CHECK(std::abs(lat - latref) < 1e-8);
    CHECK(std::abs(lon - lonref) < 1e-8);
    CHECK(std::abs(alt - altref) < 1e-8);

    CHECK(lat == doctest::Approx(latref));
    CHECK(lon == doctest::Approx(lonref));
    CHECK(alt == doctest::Approx(altref));
}

TEST_CASE("ecef2geodetic") {
    double latref = radians(30.14988205);
    double lonref = radians(91.38733072);
    double altref = 4031.0;

    double x, y, z;
    cppmap3d::geodetic2ecef(latref, lonref, altref, x, y, z);
    double lat, lon, alt;
    cppmap3d::ecef2geodetic(x, y, z, lat, lon, alt);

    CHECK(std::abs(lat - latref) < 1e-8);
    CHECK(std::abs(lon - lonref) < 1e-8);
    CHECK(std::abs(alt - altref) < 1e-8);

    CHECK(lat == doctest::Approx(latref));
    CHECK(lon == doctest::Approx(lonref));
    CHECK(alt == doctest::Approx(altref));

    cppmap3d::geodetic2ecef(latref, lonref, altref - 5000.0, x, y, z);
    cppmap3d::ecef2geodetic(x, y, z, lat, lon, alt);

    CHECK(std::abs(lat - latref) < 1e-8);
    CHECK(std::abs(lon - lonref) < 1e-8);
    CHECK(std::abs(alt - (altref - 5000.0)) < 1e-8);

    CHECK(lat == doctest::Approx(latref));
    CHECK(lon == doctest::Approx(lonref));
    CHECK(alt == doctest::Approx((altref - 5000.0)));
}

TEST_CASE("ecef2enu") {
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;
    double eref = 1.862775210000000e+02;
    double nref = 2.868422200000000e+02;
    double uref = 9.396926200000000e+02;

    double x, y, z;
    cppmap3d::enu2ecef(eref, nref, uref, lat0, lon0, alt0, x, y, z);
    double e, n, u;
    cppmap3d::ecef2enu(x, y, z, lat0, lon0, alt0, e, n, u);

    CHECK(std::abs(e - eref) < 1e-3);
    CHECK(std::abs(n - nref) < 1e-3);
    CHECK(std::abs(u - uref) < 1e-3);

    CHECK(e == doctest::Approx(eref));
    CHECK(n == doctest::Approx(nref));
    CHECK(u == doctest::Approx(uref));
}

TEST_CASE("ecef2aer") {
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;

    double azref = radians(33.0);
    double elref = radians(70.0);
    double rangeref = 1000.0;

    double x, y, z;
    cppmap3d::aer2ecef(azref, elref, rangeref, lat0, lon0, alt0, x, y, z);
    double az, el, range;
    cppmap3d::ecef2aer(x, y, z, lat0, lon0, alt0, az, el, range);

    CHECK(std::abs(az - azref) < 1e-3);
    CHECK(std::abs(el - elref) < 1e-3);
    CHECK(std::abs(range - rangeref) < 1e-3);

    CHECK(az == doctest::Approx(azref));
    CHECK(el == doctest::Approx(elref));
    CHECK(range == doctest::Approx(rangeref));
}

TEST_CASE("ned2geodetic") {
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;
    double e = 0.0;
    double n = 0.0;
    double d = 1.0;

    double latref = radians(41.999999999999993);
    double lonref = radians(-82.0);
    double altref = 1.990000000007368e+02;

    double lat, lon, alt;
    cppmap3d::ned2geodetic(n, e, d, lat0, lon0, alt0, lat, lon, alt);

    CHECK(std::abs(lat - latref) < 1e-8);
    CHECK(std::abs(lon - lonref) < 1e-8);
    CHECK(std::abs(alt - altref) < 1e-8);

    CHECK(lat == doctest::Approx(latref));
    CHECK(lon == doctest::Approx(lonref));
    CHECK(alt == doctest::Approx(altref));
}

TEST_CASE("geodetic2ned") {
    double lat = radians(41.999999999999993);
    double lon = radians(-82.0);
    double alt = 1.990000000007368e+02;
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;

    double eref = 0.0;
    double nref = 0.0;
    double dref = 1.0;

    double n, e, d;
    cppmap3d::geodetic2ned(lat, lon, alt, lat0, lon0, alt0, n, e, d);

    CHECK(std::abs(n - nref) < 1e-3);
    CHECK(std::abs(e - eref) < 1e-3);
    CHECK(std::abs(d - dref) < 1e-3);

    CHECK(e == doctest::Approx(eref));
    CHECK(n == doctest::Approx(nref));
    CHECK(d == doctest::Approx(dref));
}

TEST_CASE("aer2ned") {
    double az = radians(33.0);
    double el = radians(70.0);
    double range = 1000.0;

    double eref = 1.862775208165935e+02;
    double nref = 2.868422278517140e+02;
    double dref = -9.396926207859083e+02;

    double n, e, d;
    cppmap3d::aer2ned(az, el, range, n, e, d);

    CHECK(std::abs(n - nref) < 1e-3);
    CHECK(std::abs(e - eref) < 1e-3);
    CHECK(std::abs(d - dref) < 1e-3);

    CHECK(e == doctest::Approx(eref));
    CHECK(n == doctest::Approx(nref));
    CHECK(d == doctest::Approx(dref));
}

TEST_CASE("ned2aer") {
    double azref = radians(33.0);
    double elref = radians(70.0);
    double rangeref = 1000.0;

    double e = 1.862775208165935e+02;
    double n = 2.868422278517140e+02;
    double d = -9.396926207859083e+02;

    double az, el, range;
    cppmap3d::ned2aer(n, e, d, az, el, range);

    CHECK(std::abs(az - azref) < 1e-6);
    CHECK(std::abs(el - elref) < 1e-6);
    CHECK(std::abs(range - rangeref) < 1e-3);

    CHECK(az == doctest::Approx(azref));
    CHECK(el == doctest::Approx(elref));
    CHECK(range == doctest::Approx(rangeref));
}

TEST_CASE("ned2ecef") {
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;
    double e = 1.862775210000000e+02;
    double n = 2.868422200000000e+02;
    double d = -9.396926200000000e+02;

    double xref = 6.609301927610815e+05;
    double yref = -4.701424222957011e+06;
    double zref = 4.246579604632881e+06;

    double x, y, z;
    cppmap3d::ned2ecef(n, e, d, lat0, lon0, alt0, x, y, z);

    CHECK(std::abs(x - xref) < 1e-3);
    CHECK(std::abs(y - yref) < 1e-3);
    CHECK(std::abs(z - zref) < 1e-3);

    CHECK(x == doctest::Approx(xref));
    CHECK(y == doctest::Approx(yref));
    CHECK(z == doctest::Approx(zref));
}

TEST_CASE("ellipsoid_references") {
    cppmap3d::Ellipsoid ellipsoid = cppmap3d::Ellipsoid::WGS84;
    double a = cppmap3d::internal::getMajor(ellipsoid);
    double f = cppmap3d::internal::getFlattening(ellipsoid);
    double b = cppmap3d::internal::getMinor(ellipsoid);
    double e = cppmap3d::internal::getSquaredEccentricity(ellipsoid);

    CHECK(std::abs(a - 6378137.0) < 1E-6);
    CHECK(std::abs(b - 6356752.314245) < 1E-6);
    CHECK(std::abs(1.0 / f - 298.257223563) < 1E-6);
    CHECK(std::abs(e - 6.6943799E-3) < 1E-6);
    CHECK(a == doctest::Approx(6378137.0));
    CHECK(b == doctest::Approx(6356752.314245));
    CHECK(1.0 / f == doctest::Approx(298.257223563));
    CHECK(e == doctest::Approx(6.6943799E-3));

    ellipsoid = cppmap3d::Ellipsoid::WGS72;
    a = cppmap3d::internal::getMajor(ellipsoid);
    f = cppmap3d::internal::getFlattening(ellipsoid);
    b = cppmap3d::internal::getMinor(ellipsoid);
    e = cppmap3d::internal::getSquaredEccentricity(ellipsoid);

    CHECK(std::abs(b - a * (1.0 - f)) < 1E-6);
    CHECK((e - (f * (2.0 - f))) < 1E-6);
    CHECK(b == doctest::Approx(a * (1.0 - f)));
    CHECK(e == doctest::Approx(f * (2.0 - f)));

    ellipsoid = cppmap3d::Ellipsoid::WGS66;
    a = cppmap3d::internal::getMajor(ellipsoid);
    f = cppmap3d::internal::getFlattening(ellipsoid);
    b = cppmap3d::internal::getMinor(ellipsoid);
    e = cppmap3d::internal::getSquaredEccentricity(ellipsoid);

    CHECK(std::abs(b - a * (1.0 - f)) < 1E-6);
    CHECK((e - (f * (2.0 - f))) < 1E-6);
    CHECK(b == doctest::Approx(a * (1.0 - f)));
    CHECK(e == doctest::Approx(f * (2.0 - f)));

    CHECK(std::abs(b - a * (1.0 - f)) < 1E-6);
    CHECK((e - (f * (2.0 - f))) < 1E-6);
    CHECK(b == doctest::Approx(a * (1.0 - f)));
    CHECK(e == doctest::Approx(f * (2.0 - f)));

    ellipsoid = cppmap3d::Ellipsoid::PZ9011;
    a = cppmap3d::internal::getMajor(ellipsoid);
    f = cppmap3d::internal::getFlattening(ellipsoid);
    b = cppmap3d::internal::getMinor(ellipsoid);
    e = cppmap3d::internal::getSquaredEccentricity(ellipsoid);

    CHECK(std::abs(a - 6378136.0) < 1E-6);
    CHECK(std::abs(b - a * (1.0 - f)) < 1E-6);
    CHECK(std::abs(1.0 / f - 298.257839303) < 1E-6);
    CHECK((e - (f * (2.0 - f))) < 1E-6);
    CHECK(a == doctest::Approx(6378136.0));
    CHECK(b == doctest::Approx(a * (1.0 - f)));
    CHECK(1.0 / f == doctest::Approx(298.257839303));
    CHECK(e == doctest::Approx(f * (2.0 - f)));

    ellipsoid = cppmap3d::Ellipsoid::GRS80;
    a = cppmap3d::internal::getMajor(ellipsoid);
    f = cppmap3d::internal::getFlattening(ellipsoid);
    b = cppmap3d::internal::getMinor(ellipsoid);
    e = cppmap3d::internal::getSquaredEccentricity(ellipsoid);

    CHECK(std::abs(a - 6378137.0) < 1E-6);
    CHECK(std::abs(b - a * (1.0 - f)) < 1E-6);
    CHECK(std::abs(1.0 / f - 298.257222101) < 1E-6);
    CHECK((e - (f * (2.0 - f))) < 1E-6);
    CHECK(a == doctest::Approx(6378137.0));
    CHECK(b == doctest::Approx(a * (1.0 - f)));
    CHECK(1.0 / f == doctest::Approx(298.257222101));
    CHECK(e == doctest::Approx(f * (2.0 - f)));

    ellipsoid = cppmap3d::Ellipsoid::Bessel;
    a = cppmap3d::internal::getMajor(ellipsoid);
    f = cppmap3d::internal::getFlattening(ellipsoid);
    b = cppmap3d::internal::getMinor(ellipsoid);
    e = cppmap3d::internal::getSquaredEccentricity(ellipsoid);

    CHECK(std::abs(b - a * (1.0 - f)) < 1E-6);
    CHECK((e - (f * (2.0 - f))) < 1E-6);
    CHECK(b == doctest::Approx(a * (1.0 - f)));
    CHECK(e == doctest::Approx(f * (2.0 - f)));

    ellipsoid = cppmap3d::Ellipsoid::International1967;
    a = cppmap3d::internal::getMajor(ellipsoid);
    f = cppmap3d::internal::getFlattening(ellipsoid);
    b = cppmap3d::internal::getMinor(ellipsoid);
    e = cppmap3d::internal::getSquaredEccentricity(ellipsoid);

    CHECK(std::abs(b - a * (1.0 - f)) < 1E-6);
    CHECK((e - (f * (2.0 - f))) < 1E-6);
    CHECK(b == doctest::Approx(a * (1.0 - f)));
    CHECK(e == doctest::Approx(f * (2.0 - f)));

    ellipsoid = cppmap3d::Ellipsoid::Airy;
    a = cppmap3d::internal::getMajor(ellipsoid);
    f = cppmap3d::internal::getFlattening(ellipsoid);
    b = cppmap3d::internal::getMinor(ellipsoid);
    e = cppmap3d::internal::getSquaredEccentricity(ellipsoid);

    CHECK(std::abs(b - a * (1.0 - f)) < 1E-6);
    CHECK((e - (f * (2.0 - f))) < 1E-6);
    CHECK(b == doctest::Approx(a * (1.0 - f)));
    CHECK(e == doctest::Approx(f * (2.0 - f)));
}

TEST_CASE("ecef2ned") {
    double lat0 = radians(42.0);
    double lon0 = radians(-82.0);
    double alt0 = 200.0;
    double eref = 1.862775210000000e+02;
    double nref = 2.868422200000000e+02;
    double dref = -9.396926200000000e+02;

    double x, y, z;
    double n, e, d;
    cppmap3d::ned2ecef(nref, eref, dref, lat0, lon0, alt0, x, y, z);
    cppmap3d::ecef2ned(x, y, z, lat0, lon0, alt0, n, e, d);

    CHECK(std::abs(n - nref) < 1e-3);
    CHECK(std::abs(e - eref) < 1e-3);
    CHECK(std::abs(d - dref) < 1e-3);

    CHECK(e == doctest::Approx(eref));
    CHECK(n == doctest::Approx(nref));
    CHECK(d == doctest::Approx(dref));
}
