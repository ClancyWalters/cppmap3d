#ifndef CPPMAP3D_HEADER_GUARD
#define CPPMAP3D_HEADER_GUARD

#include <cmath>
#include <stdexcept>

/**
 * C++ 3D geographic coordinate conversions
 */
namespace cppmap3d {

/**
 * @brief Ellipsoids used to approximate the shape of the earth.
 *
 * @see
 * [pymap3d ellipsoid
 * source](https://github.com/geospace-code/pymap3d/blob/main/src/pymap3d/ellipsoid.py)
 */
enum class Ellipsoid {
    Maupertuis,
    Plessis,
    Everest1830,
    Everest1830Modified,
    Everest1967,
    Airy,
    Bessel,
    Clarke1866,
    Clarke1878,
    Clarke1860,
    Helmert,
    Hayford,
    International1924,
    Krassovsky1940,
    WGS66,
    Australian,
    International1967,
    GRS67,
    SA1969,
    WGS72,
    GRS80,
    WGS84,
    WGS84Mean,
    IERS1989,
    PZ9011,
    IERS2003,
    GSK2011,
    Mercury,
    Venus,
    Moon,
    Mars,
    Jupyter,
    Io,
    Saturn,
    Uranus,
    Neptune,
    Pluto
};

namespace internal {
/**
 * @brief Returns the semi-major axis length of the specified ellipsoid.
 * @param ellipsoid The ellipsoid for which the semi-major axis length is
 * requested.
 * @return The semi-major axis length of the specified ellipsoid in meters.
 * @note The semi-major axis represents half of the longest diameter of the
 * ellipsoid.
 */
inline constexpr double getMajor(Ellipsoid ellipsoid) {
    switch (ellipsoid) {
        case cppmap3d::Ellipsoid::Maupertuis:
            return 6397300.0;
        case cppmap3d::Ellipsoid::Plessis:
            return 6376523.0;
        case cppmap3d::Ellipsoid::Everest1830:
            return 6377299.365;
        case cppmap3d::Ellipsoid::Everest1830Modified:
            return 6377304.063;
        case cppmap3d::Ellipsoid::Everest1967:
            return 6377298.556;
        case cppmap3d::Ellipsoid::Airy:
            return 6377563.396;
        case cppmap3d::Ellipsoid::Bessel:
            return 6377397.155;
        case cppmap3d::Ellipsoid::Clarke1866:
            return 6378206.4;
        case cppmap3d::Ellipsoid::Clarke1878:
            return 6378190.0;
        case cppmap3d::Ellipsoid::Clarke1860:
            return 6378249.145;
        case cppmap3d::Ellipsoid::Helmert:
            return 6378200.0;
        case cppmap3d::Ellipsoid::Hayford:
            return 6378388.0;
        case cppmap3d::Ellipsoid::International1924:
            return 6378388.0;
        case cppmap3d::Ellipsoid::Krassovsky1940:
            return 6378245.0;
        case cppmap3d::Ellipsoid::WGS66:
            return 6378145.0;
        case cppmap3d::Ellipsoid::Australian:
            return 6378160.0;
        case cppmap3d::Ellipsoid::International1967:
            return 6378157.5;
        case cppmap3d::Ellipsoid::GRS67:
            return 6378160.0;
        case cppmap3d::Ellipsoid::SA1969:
            return 6378160.0;
        case cppmap3d::Ellipsoid::WGS72:
            return 6378135.0;
        case cppmap3d::Ellipsoid::GRS80:
            return 6378137.0;
        case cppmap3d::Ellipsoid::WGS84:
            return 6378137.0;
        case cppmap3d::Ellipsoid::WGS84Mean:
            return 6371008.7714;
        case cppmap3d::Ellipsoid::IERS1989:
            return 6378136.0;
        case cppmap3d::Ellipsoid::PZ9011:
            return 6378136.0;
        case cppmap3d::Ellipsoid::IERS2003:
            return 6378136.6;
        case cppmap3d::Ellipsoid::GSK2011:
            return 6378136.5;
        case cppmap3d::Ellipsoid::Mercury:
            return 2440500.0;
        case cppmap3d::Ellipsoid::Venus:
            return 6051800.0;
        case cppmap3d::Ellipsoid::Moon:
            return 1738100.0;
        case cppmap3d::Ellipsoid::Mars:
            return 3396900.0;
        case cppmap3d::Ellipsoid::Jupyter:
            return 71492000.0;
        case cppmap3d::Ellipsoid::Io:
            return 1829.7;
        case cppmap3d::Ellipsoid::Saturn:
            return 60268000.0;
        case cppmap3d::Ellipsoid::Uranus:
            return 25559000.0;
        case cppmap3d::Ellipsoid::Neptune:
            return 24764000.0;
        case cppmap3d::Ellipsoid::Pluto:
            return 1188000.0;
        default:
            return getMajor(cppmap3d::Ellipsoid::WGS84);
    }
}

/**
 * @brief Returns the semi-minor axis length of the specified ellipsoid.
 * @param ellipsoid The ellipsoid for which the semi-minor axis length is
 * requested.
 * @return The semi-minor axis length of the specified ellipsoid in meters.
 * @note The semi-minor axis represents half of the shortest diameter of the
 * ellipsoid.
 */
inline constexpr double getMinor(Ellipsoid ellipsoid) {
    switch (ellipsoid) {
        case cppmap3d::Ellipsoid::Maupertuis:
            return 6363806.283;
        case cppmap3d::Ellipsoid::Plessis:
            return 6355862.9333;
        case cppmap3d::Ellipsoid::Everest1830:
            return 6356098.359;
        case cppmap3d::Ellipsoid::Everest1830Modified:
            return 6356103.039;
        case cppmap3d::Ellipsoid::Everest1967:
            return 6356097.55;
        case cppmap3d::Ellipsoid::Airy:
            return 6356256.909;
        case cppmap3d::Ellipsoid::Bessel:
            return 6356078.963;
        case cppmap3d::Ellipsoid::Clarke1866:
            return 6356583.8;
        case cppmap3d::Ellipsoid::Clarke1878:
            return 6356456.0;
        case cppmap3d::Ellipsoid::Clarke1860:
            return 6356514.87;
        case cppmap3d::Ellipsoid::Helmert:
            return 6356818.17;
        case cppmap3d::Ellipsoid::Hayford:
            return 6356911.946;
        case cppmap3d::Ellipsoid::International1924:
            return 6356911.946;
        case cppmap3d::Ellipsoid::Krassovsky1940:
            return 6356863.019;
        case cppmap3d::Ellipsoid::WGS66:
            return 6356759.769;
        case cppmap3d::Ellipsoid::Australian:
            return 6356774.719;
        case cppmap3d::Ellipsoid::International1967:
            return 6356772.2;
        case cppmap3d::Ellipsoid::GRS67:
            return 6356774.516;
        case cppmap3d::Ellipsoid::SA1969:
            return 6356774.719;
        case cppmap3d::Ellipsoid::WGS72:
            return 6356750.52001609;
        case cppmap3d::Ellipsoid::GRS80:
            return 6356752.31414036;
        case cppmap3d::Ellipsoid::WGS84:
            return 6356752.31424518;
        case cppmap3d::Ellipsoid::WGS84Mean:
            return 6371008.7714;
        case cppmap3d::Ellipsoid::IERS1989:
            return 6356751.302;
        case cppmap3d::Ellipsoid::PZ9011:
            return 6356751.3618;
        case cppmap3d::Ellipsoid::IERS2003:
            return 6356751.9;
        case cppmap3d::Ellipsoid::GSK2011:
            return 6356751.758;
        case cppmap3d::Ellipsoid::Mercury:
            return 2438300.0;
        case cppmap3d::Ellipsoid::Venus:
            return 6051800.0;
        case cppmap3d::Ellipsoid::Moon:
            return 1736000.0;
        case cppmap3d::Ellipsoid::Mars:
            return 3376097.80585952;
        case cppmap3d::Ellipsoid::Jupyter:
            return 66770054.3475922;
        case cppmap3d::Ellipsoid::Io:
            return 1815.8;
        case cppmap3d::Ellipsoid::Saturn:
            return 54364301.5271271;
        case cppmap3d::Ellipsoid::Uranus:
            return 24973000.0;
        case cppmap3d::Ellipsoid::Neptune:
            return 24341000.0;
        case cppmap3d::Ellipsoid::Pluto:
            return 1188000.0;
        default:
            return getMinor(cppmap3d::Ellipsoid::WGS84);
    }
}

/**
 * @brief Returns the flattening factor of the specified ellipsoid.
 * @param ellipsoid The ellipsoid for which the flattening factor is requested.
 * @return The flattening factor, a dimensionless quantity representing the
 * ellipsoid's deviation from a perfect sphere.
 * @note Flattening is calculated as (a - b) / a, where 'a' is the semi-major
 * axis and 'b' is the semi-minor axis.
 */
inline constexpr double getFlattening(Ellipsoid ellipsoid) {
    double major = getMajor(ellipsoid);
    double minor = getMinor(ellipsoid);
    return (major - minor) / major;
}

/**
 * @brief Returns the squared eccentricity of the specified ellipsoid.
 * @param ellipsoid The ellipsoid for which the squared eccentricity is
 * requested.
 * @return The squared eccentricity, a dimensionless quantity expressing the
 * ellipsoid's deviation from a perfect sphere.
 * @note Squared eccentricity is calculated as (a^2 - b^2) / a^2, where 'a' is
 * the semi-major axis and 'b' is the semi-minor axis.
 */
inline constexpr double getSquaredEccentricity(Ellipsoid ellipsoid) {
    double major = getMajor(ellipsoid);
    double minor = getMinor(ellipsoid);
    return ((major * major) - (minor * minor)) / (major * major);
}

/**
 * @brief Internal function used in conversions certain coordinate conversions
 */
inline void enu2uvw(
    double et,
    double nt,
    double up,
    double lat,
    double lon,
    double& out_u,
    double& out_v,
    double& out_w
) {
    double t = std::cos(lat) * up - std::sin(lat) * nt;
    out_u = std::cos(lon) * t - std::sin(lon) * et;
    out_v = std::sin(lon) * t + std::cos(lon) * et;
    out_w = std::sin(lat) * up + std::cos(lat) * nt;
}

/**
 * @brief Internal function used in conversions certain coordinate conversions
 */
inline void uvw2enu(
    double u,
    double v,
    double w,
    double lat,
    double lon,
    double& out_east,
    double& out_north,
    double& out_up
) {
    double t = std::cos(lon) * u + std::sin(lon) * v;
    out_east = -1 * std::sin(lon) * u + std::cos(lon) * v;
    out_north = -1 * std::sin(lat) * t + std::cos(lat) * w;
    out_up = std::cos(lat) * t + std::sin(lat) * w;
}

/**
 * @brief Converts Earth-Centered, Earth-Fixed (ECEF) coordinates to geodetic
 * coordinates.
 *
 * @param x The ECEF x-coordinate, in meters.
 * @param y The ECEF y-coordinate, in meters.
 * @param z The ECEF z-coordinate, in meters.
 * @param[out] out_lat The latitude, in radians (output parameter).
 * @param[out] out_lon The longitude, in radians (output parameter).
 * @param[out] out_alt The altitude, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 * @note This implementation is ported from pymap3d and is based on: You,
 * Rey-Jer. (2000). Transformation of Cartesian to Geodetic Coordinates without
 * Iterations. Journal of Surveying Engineering. doi: 10.1061/(ASCE)0733-9453.
 *
 * @note The latitude and longitude are returned in radians, and the altitude
 * is provided in meters.
 *
 */
inline void ecef2geodetic_you(
    double x,
    double y,
    double z,
    double& out_lat,
    double& out_lon,
    double& out_alt,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double major = internal::getMajor(ellipsoid);
    double minor = internal::getMinor(ellipsoid);

    double r = std::sqrt(x * x + y * y + z * z);
    double e = std::sqrt(major * major - minor * minor);
    double var = r * r - e * e;
    double u =
        std::sqrt(0.5 * var + 0.5 * std::sqrt(var * var + 4.0 * e * e * z * z));

    double q = std::sqrt(x * x + y * y);
    double hu_e = std::sqrt(u * u + e * e);

    // it's possible for this to be nan since u
    double beta = std::atan(hu_e / u * z / q);

    if (std::isnan(beta)) {
        if (std::abs(z) < 1.0e-9) {
            beta = 0;
        } else if (z > 0) {
            beta = 3.14159265358979311599796346854 / 2;
        } else {
            beta = -3.14159265358979311599796346854 / 2;
        }
    }

    double eps = ((minor * u - major * hu_e + e * e) * std::sin(beta)) /
                 (major * hu_e / std::cos(beta) - e * e * std::cos(beta));

    beta += eps;

    out_lat = std::atan(major / minor * std::tan(beta));
    out_lon = std::atan2(y, x);

    double v1 = z - minor * std::sin(beta);
    double v2 = q - major * std::cos(beta);

    if ((x * x / major / major) + (y * y / major / major) +
            (z * z / minor / minor) <
        1.0) {
        out_alt = -1 * std::sqrt(v1 * v1 + v2 * v2);
    } else {
        out_alt = std::sqrt(v1 * v1 + v2 * v2);
    }
}

/**
 * @brief Converts Earth-Centered, Earth-Fixed (ECEF) coordinates to geodetic
 * coordinates using the Olson (1996) algorithm.
 *
 * @param x The ECEF x-coordinate, in meters.
 * @param y The ECEF y-coordinate, in meters.
 * @param z The ECEF z-coordinate, in meters.
 * @param[out] out_lat The latitude, in radians (output parameter).
 * @param[out] out_lon The longitude, in radians (output parameter).
 * @param[out] out_alt The altitude, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 * @note This function implements the Olson (1996) algorithm, sourced from
 * [this algorithm
 * comparison](https://github.com/planet36/ecef-geodetic/blob/main/olson_1996/olson_1996.c)
 *
 * @throws std::domain_error if the ECEF coordinates are close to the center of
 * the Earth.
 */
inline void ecef2geodetic_olson(
    double x,
    double y,
    double z,
    double& out_lat,
    double& out_lon,
    double& out_alt,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double zp, w2, w, z2, r2, r, s2, c2, s, c, ss;
    double g, rg, rf, u, v, m, f, p;
    double a = internal::getMajor(ellipsoid);
    double e2 = internal::getSquaredEccentricity(ellipsoid);
    double a1 = a * e2;
    double a2 = a1 * a1;
    double a3 = a1 * e2 / 2.0;
    double a4 = (5.0 / 2.0) * a2;
    double a5 = a1 + a3;
    double a6 = 1 - e2;

    zp = std::fabs(z);
    w2 = x * x + y * y;
    w = std::sqrt(w2);
    z2 = z * z;
    r2 = w2 + z2;
    r = std::sqrt(r2);
    if (r < 100000.) {
        out_lat = 0.;
        out_lon = 0.;
        out_alt = -1.e7;
        throw std::domain_error(
            "Cannot calculate geodetic of ecef close to center of earth"
        );
    }
    out_lon = std::atan2(y, x);
    s2 = z2 / r2;
    c2 = w2 / r2;
    u = a2 / r;
    v = a3 - a4 / r;
    if (c2 > .3) {
        s = (zp / r) * (1. + c2 * (a1 + u + s2 * v) / r);
        out_lat = std::asin(s);
        ss = s * s;
        c = std::sqrt(1. - ss);
    } else {
        c = (w / r) * (1. - s2 * (a5 - u - c2 * v) / r);
        out_lat = std::acos(c);
        ss = 1. - c * c;
        s = std::sqrt(ss);
    }
    g = 1. - e2 * ss;
    rg = a / std::sqrt(g);
    rf = a6 * rg;
    u = w - rg * c;
    v = zp - rf * s;
    f = c * u + s * v;
    m = c * v - s * u;
    p = m / (rf / g + f);
    out_lat = out_lat + p;
    out_alt = f + m * p / 2.;
    if (z < 0.) {
        out_lat = -out_lat;
    }
}

}  // namespace internal

/**
 * @brief Converts Earth-Centered, Earth-Fixed (ECEF) coordinates to geodetic
 * coordinates.
 *
 * @param x The ECEF x-coordinate, in meters.
 * @param y The ECEF y-coordinate, in meters.
 * @param z The ECEF z-coordinate, in meters.
 * @param[out] out_lat The latitude, in radians (output parameter).
 * @param[out] out_lon The longitude, in radians (output parameter).
 * @param[out] out_alt The altitude, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 * @note Defaults to pymap3d's algorithm from you (2000).
 * Olson (1996) algorithm is selected at compile time by defining
 * the 'CPPMAP3D_ECEF2GEODETIC_OLSON' macro.
 *
 * @throws std::domain_error if the ECEF coordinates are close to the center of
 * the Earth and CPPMAP3D_ECEF2GEODETIC_OLSON is defined.
 */
inline void ecef2geodetic(
    double x,
    double y,
    double z,
    double& out_lat,
    double& out_lon,
    double& out_alt,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
#ifndef CPPMAP3D_ECEF2GEODETIC_OLSON
    internal::ecef2geodetic_you(x, y, z, out_lat, out_lon, out_alt, ellipsoid);
#endif
#ifdef CPPMAP3D_ECEF2GEODETIC_OLSON
    internal::ecef2geodetic_olson(
        x,
        y,
        z,
        out_lat,
        out_lon,
        out_alt,
        ellipsoid
    );
#endif
}

/**
 * @brief Converts Azimuth, Elevation, and Range (AER) coordinates to East,
 * North, Up (ENU) coordinates.
 *
 * @param az The ECEF azimuth angle, in radians.
 * @param el The elevation angle, in radians.
 * @param range The range distance, in meters.
 * @param[out] out_e The East coordinate, in meters (output parameter).
 * @param[out] out_n The North coordinate, in meters (output parameter).
 * @param[out] out_u The Up coordinate, in meters (output parameter).
 *
 * @throws std::domain_error if range is negative
 */
inline void aer2enu(
    double az,
    double el,
    double range,
    double& out_e,
    double& out_n,
    double& out_u
) {
    if (range < 0) {
        throw std::domain_error("range should not be negative");
    }

    auto r = range * std::cos(el);
    out_e = r * std::sin(az);
    out_n = r * std::cos(az);
    out_u = range * std::sin(el);
}

/**
 * @brief Converts East, North, Up (ENU) coordinates to Azimuth, Elevation, and
 * Range (AER) coordinates.
 *
 * @param east The East coordinate, in meters.
 * @param north The North coordinate, in meters.
 * @param up The Up coordinate, in meters.
 * @param[out] out_az The ECEF azimuth angle, in radians (output parameter).
 * @param[out] out_el The elevation angle, in radians (output parameter).
 * @param[out] out_range The range distance, in meters (output parameter).
 *
 */
inline void enu2aer(
    double east,
    double north,
    double up,
    double& out_az,
    double& out_el,
    double& out_range
) {
    // 1 millimeter precision for singularity stability - see pymap3d PR#42
    if (std::abs(east) < 1.0e-3) {
        east = 0;
    }
    if (std::abs(north) < 1.0e-3) {
        north = 0;
    }
    if (std::abs(up) < 1.0e-3) {
        up = 0;
    }

    double r = std::hypot(east, north);

    out_range = std::hypot(r, up);
    out_el = std::atan2(up, r);
    const auto pi2 = (2 * 3.14159265358979311599796346854);

    out_az = std::fmod(pi2 + std::fmod(std::atan2(east, north), pi2), pi2);
}

/**
 * @brief Converts Azimuth, Elevation, and Range (AER) coordinates to North,
 * East, Down (NED) coordinates.
 *
 * @param az The ECEF azimuth angle, in radians.
 * @param el The elevation angle, in radians.
 * @param range The range distance, in meters.
 * @param[out] out_north The North coordinate, in meters (output parameter).
 * @param[out] out_east The East coordinate, in meters (output parameter).
 * @param[out] out_down The Down coordinate, in meters (output parameter).
 *
 */
inline void aer2ned(
    double az,
    double el,
    double range,
    double& out_north,
    double& out_east,
    double& out_down
) {
    aer2enu(az, el, range, out_east, out_north, out_down);
    out_down = out_down * -1;
}

/**
 * @brief Converts geodetic coordinates to Earth-Centered, Earth-Fixed (ECEF)
 * coordinates.
 *
 * @param lat The latitude, in radians.
 * @param lon The longitude, in radians.
 * @param alt The altitude, in meters.
 * @param[out] out_x The ECEF x-coordinate, in meters (output
 * parameter).
 * @param[out] out_y The ECEF y-coordinate, in meters (output
 * parameter).
 * @param[out] out_z The ECEF z-coordinate, in meters (output
 * parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 * @throws std::domain_error if latitude is not within -pi/2 and pi/2 inclusive
 */
inline void geodetic2ecef(
    double lat,
    double lon,
    double alt,
    double& out_x,
    double& out_y,
    double& out_z,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    const double pi = 3.14159265358979311599796346854;
    if (std::abs(lat) > pi / 2) {
        throw std::domain_error("-pi/2 <= latitude <= pi/2");
    }

    double major = internal::getMajor(ellipsoid);
    double minor = internal::getMinor(ellipsoid);
    double se = internal::getSquaredEccentricity(ellipsoid);

    double n = major / std::sqrt(((1.0 - se * std::sin(lat) * std::sin(lat))));

    out_x = (n + alt) * std::cos(lat) * std::cos(lon);
    out_y = (n + alt) * std::cos(lat) * std::sin(lon);
    out_z = (n * (minor / major) * (minor / major) + alt) * std::sin(lat);
}

/**
 * @brief Converts Earth-Centered, Earth-Fixed (ECEF) coordinates to East,
 * North, Up (ENU) coordinates.
 *
 * @param x The ECEF x-coordinate, in meters.
 * @param y The ECEF y-coordinate, in meters.
 * @param z The ECEF z-coordinate, in meters.
 * @param lat The latitude, in radians.
 * @param lon The longitude, in radians.
 * @param alt The altitude, in meters.
 * @param[out] out_east The East coordinate, in meters (output parameter).
 * @param[out] out_north The North coordinate, in meters (output parameter).
 * @param[out] out_up The Up coordinate, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void ecef2enu(
    double x,
    double y,
    double z,
    double lat,
    double lon,
    double alt,
    double& out_east,
    double& out_north,
    double& out_up,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double x0, y0, z0;
    geodetic2ecef(lat, lon, alt, x0, y0, z0, ellipsoid);
    internal::uvw2enu(
        x - x0,
        y - y0,
        z - z0,
        lat,
        lon,
        out_east,
        out_north,
        out_up
    );
}

/**
 * @brief Converts Earth-Centered, Earth-Fixed (ECEF) coordinates to Azimuth,
 * Elevation, and Range (AER) coordinates.
 *
 * @param x The ECEF x-coordinate, in meters.
 * @param y The ECEF y-coordinate, in meters.
 * @param z The ECEF z-coordinate, in meters.
 * @param lat The latitude, in radians.
 * @param lon The longitude, in radians.
 * @param alt The altitude, in meters.
 * @param[out] out_az The ECEF azimuth angle, in radians (output parameter).
 * @param[out] out_el The elevation angle, in radians (output parameter).
 * @param[out] out_range The range distance, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void ecef2aer(
    double x,
    double y,
    double z,
    double lat,
    double lon,
    double alt,
    double& out_az,
    double& out_el,
    double& out_range,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double e, n, u;
    ecef2enu(x, y, z, lat, lon, alt, e, n, u, ellipsoid);
    enu2aer(e, n, u, out_az, out_el, out_range);
}

/**
 * @brief Converts Earth-Centered, Earth-Fixed (ECEF) coordinates to North,
 * East, Down (NED) coordinates.
 *
 * @param x The ECEF x-coordinate, in meters.
 * @param y The ECEF y-coordinate, in meters.
 * @param z The ECEF z-coordinate, in meters.
 * @param lat The latitude, in radians.
 * @param lon The longitude, in radians.
 * @param alt The altitude, in meters.
 * @param[out] out_north The North coordinate, in meters (output parameter).
 * @param[out] out_east The East coordinate, in meters (output parameter).
 * @param[out] out_down The Down coordinate, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void ecef2ned(
    double x,
    double y,
    double z,
    double lat,
    double lon,
    double alt,
    double& out_north,
    double& out_east,
    double& out_down,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    ecef2enu(x, y, z, lat, lon, alt, out_east, out_north, out_down, ellipsoid);
    out_down = out_down * -1;
}

/**
 * @brief Converts East, North, Up (ENU) coordinates to Earth-Centered,
 * Earth-Fixed (ECEF) coordinates.
 *
 * @param east The East coordinate, in meters.
 * @param north The North coordinate, in meters.
 * @param up The Up coordinate, in meters.
 * @param lat The latitude, in radians.
 * @param lon The longitude, in radians.
 * @param alt The altitude, in meters.
 * @param[out] out_x The ECEF x-coordinate, in meters (output parameter).
 * @param[out] out_y The ECEF y-coordinate, in meters (output parameter).
 * @param[out] out_z The ECEF z-coordinate, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void enu2ecef(
    double east,
    double north,
    double up,
    double lat,
    double lon,
    double alt,
    double& out_x,
    double& out_y,
    double& out_z,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double x, y, z;
    geodetic2ecef(lat, lon, alt, x, y, z, ellipsoid);
    internal::enu2uvw(east, north, up, lat, lon, out_x, out_y, out_z);
    out_x = out_x + x;
    out_y = out_y + y;
    out_z = out_z + z;
}

/**
 * @brief Converts East, North, Up (ENU) coordinates to geodetic coordinates.
 *
 * @param east The East coordinate, in meters.
 * @param north The North coordinate, in meters.
 * @param up The Up coordinate, in meters.
 * @param lat The latitude, in radians.
 * @param lon The longitude, in radians.
 * @param alt The altitude, in meters.
 * @param[out] out_lat The latitude, in radians (output parameter).
 * @param[out] out_lon The longitude, in radians (output parameter).
 * @param[out] out_alt The altitude, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void enu2geodetic(
    double east,
    double north,
    double up,
    double lat,
    double lon,
    double alt,
    double& out_lat,
    double& out_lon,
    double& out_alt,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double x, y, z;
    enu2ecef(east, north, up, lat, lon, alt, x, y, z, ellipsoid);
    ecef2geodetic(x, y, z, out_lat, out_lon, out_alt, ellipsoid);
}

/**
 * @brief Converts geodetic coordinates to East, North, Up (ENU) coordinates
 * relative to a reference point.
 *
 * @param lat The latitude in radians of the target.
 * @param lon The longitude in radians of the target.
 * @param alt The altitude in meters of the target.
 * @param lat0 The latitude in radians of the observer.
 * @param lon0 The longitude in radians of the observer.
 * @param alt0 The altitude in meters of the observer.
 * @param[out] out_east The East coordinate, in meters (output parameter).
 * @param[out] out_north The North coordinate, in meters (output parameter).
 * @param[out] out_up The Up coordinate, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void geodetic2enu(
    double lat,
    double lon,
    double alt,
    double lat0,
    double lon0,
    double alt0,
    double& out_east,
    double& out_north,
    double& out_up,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double x1, y1, z1;
    geodetic2ecef(lat, lon, alt, x1, y1, z1, ellipsoid);

    double x2, y2, z2;
    geodetic2ecef(lat0, lon0, alt0, x2, y2, z2, ellipsoid);

    internal::uvw2enu(
        x1 - x2,
        y1 - y2,
        z1 - z2,
        lat0,
        lon0,
        out_east,
        out_north,
        out_up
    );
}

/**
 * @brief Converts geodetic coordinates to Azimuth, Elevation, and Range (AER)
 * coordinates relative to a reference point.
 *
 * @param lat The latitude in radians of the target.
 * @param lon The longitude in radians of the target.
 * @param alt The altitude in meters of the target.
 * @param lat0 The latitude in radians of the observer.
 * @param lon0 The longitude in radians of the observer.
 * @param alt0 The altitude in meters of the observer.
 * @param[out] out_az The ECEF azimuth angle, in radians (output parameter).
 * @param[out] out_el The elevation angle, in radians (output parameter).
 * @param[out] out_range The range distance, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void geodetic2aer(
    double lat,
    double lon,
    double alt,
    double lat0,
    double lon0,
    double alt0,
    double& out_az,
    double& out_el,
    double& out_range,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double e, n, u;
    geodetic2enu(lat, lon, alt, lat0, lon0, alt0, e, n, u, ellipsoid);
    enu2aer(e, n, u, out_az, out_el, out_range);
}

/**
 * @brief Converts geodetic coordinates to North, East, Down (NED) coordinates
 * relative to a reference point.
 *
 * @param lat The latitude in radians of the target.
 * @param lon The longitude in radians of the target.
 * @param alt The altitude in meters of the target.
 * @param lat0 The latitude in radians of the observer.
 * @param lon0 The longitude in radians of the observer.
 * @param alt0 The altitude in meters of the observer.
 * @param[out] out_north The North coordinate, in meters (output parameter).
 * @param[out] out_east The East coordinate, in meters (output parameter).
 * @param[out] out_down The Down coordinate, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void geodetic2ned(
    double lat,
    double lon,
    double alt,
    double lat0,
    double lon0,
    double alt0,
    double& out_north,
    double& out_east,
    double& out_down,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    geodetic2enu(
        lat,
        lon,
        alt,
        lat0,
        lon0,
        alt0,
        out_east,
        out_north,
        out_down,
        ellipsoid
    );
    out_down = out_down * -1;
}

/**
 * @brief Converts Azimuth, Elevation, and Range (AER) coordinates to
 * Earth-Centered, Earth-Fixed (ECEF) coordinates.
 *
 * @param az The ECEF azimuth angle, in radians.
 * @param el The elevation angle, in radians.
 * @param range The range distance, in meters.
 * @param lat The latitude, in radians.
 * @param lon The longitude, in radians.
 * @param alt The altitude, in meters.
 * @param[out] out_x The ECEF x-coordinate, in meters (output
 * parameter).
 * @param[out] out_y The ECEF y-coordinate, in meters (output
 * parameter).
 * @param[out] out_z The ECEF z-coordinate, in meters (output
 * parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void aer2ecef(
    double az,
    double el,
    double range,
    double lat,
    double lon,
    double alt,
    double& out_x,
    double& out_y,
    double& out_z,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double x, y, z;
    geodetic2ecef(lat, lon, alt, x, y, z, ellipsoid);
    double e, n, u;
    aer2enu(az, el, range, e, n, u);
    double dx, dy, dz;
    internal::enu2uvw(e, n, u, lat, lon, dx, dy, dz);
    out_x = x + dx;
    out_y = y + dy;
    out_z = z + dz;
}

/**
 * @brief Converts Azimuth, Elevation, and Range (AER) coordinates to geodetic
 * coordinates.
 *
 * @param az The ECEF azimuth angle, in radians.
 * @param el The elevation angle, in radians.
 * @param range The range distance, in meters.
 * @param lat The latitude, in radians.
 * @param lon The longitude, in radians.
 * @param alt The altitude, in meters.
 * @param[out] out_lat The latitude, in radians (output parameter).
 * @param[out] out_lon The longitude, in radians (output parameter).
 * @param[out] out_alt The altitude, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void aer2geodetic(
    double az,
    double el,
    double range,
    double lat,
    double lon,
    double alt,
    double& out_lat,
    double& out_lon,
    double& out_alt,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double x, y, z;
    aer2ecef(az, el, range, lat, lon, alt, x, y, z, ellipsoid);
    ecef2geodetic(x, y, z, out_lat, out_lon, out_alt, ellipsoid);
}

/**
 * @brief Converts North, East, Down (NED) coordinates to Earth-Centered,
 * Earth-Fixed (ECEF) coordinates.
 *
 * @param north The North coordinate, in meters.
 * @param east The East coordinate, in meters.
 * @param down The Down coordinate, in meters.
 * @param lat The latitude, in radians.
 * @param lon The longitude, in radians.
 * @param alt The altitude, in meters.
 * @param[out] out_x The ECEF x-coordinate, in meters (output parameter).
 * @param[out] out_y The ECEF y-coordinate, in meters (output parameter).
 * @param[out] out_z The ECEF z-coordinate, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void ned2ecef(
    double north,
    double east,
    double down,
    double lat,
    double lon,
    double alt,
    double& out_x,
    double& out_y,
    double& out_z,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    enu2ecef(
        east,
        north,
        -1 * down,
        lat,
        lon,
        alt,
        out_x,
        out_y,
        out_z,
        ellipsoid
    );
}

/**
 * @brief Converts North, East, Down (NED) coordinates to geodetic coordinates.
 *
 * @param north The North coordinate, in meters.
 * @param east The East coordinate, in meters.
 * @param down The Down coordinate, in meters.
 * @param lat The latitude, in radians.
 * @param lon The longitude, in radians.
 * @param alt The altitude, in meters.
 * @param[out] out_lat The latitude, in radians (output parameter).
 * @param[out] out_lon The longitude, in radians (output parameter).
 * @param[out] out_alt The altitude, in meters (output parameter).
 * @param ellipsoid The ellipsoid model used for the conversion (default is
 * WGS84).
 *
 */
inline void ned2geodetic(
    double north,
    double east,
    double down,
    double lat,
    double lon,
    double alt,
    double& out_lat,
    double& out_lon,
    double& out_alt,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    enu2geodetic(
        east,
        north,
        -1 * down,
        lat,
        lon,
        alt,
        out_lat,
        out_lon,
        out_alt,
        ellipsoid
    );
}

/**
 * @brief Converts North, East, Down (NED) coordinates to Azimuth, Elevation,
 * and Range (AER) coordinates.
 *
 * @param north The North coordinate, in meters.
 * @param east The East coordinate, in meters.
 * @param down The Down coordinate, in meters.
 * @param[out] out_az The azimuth angle, in radians (output parameter).
 * @param[out] out_el The elevation angle, in radians (output parameter).
 * @param[out] out_range The range distance, in meters (output parameter).
 *
 */
inline void ned2aer(
    double north,
    double east,
    double down,
    double& out_az,
    double& out_el,
    double& out_range
) {
    enu2aer(east, north, -1 * down, out_az, out_el, out_range);
}
}  // namespace cppmap3d

#endif