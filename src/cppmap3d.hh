#ifndef CPPMAP3D_HEADER_GUARD
#define CPPMAP3D_HEADER_GUARD

#include <cmath>

namespace cppmap3d {

enum class Ellipsoid {
    /// WGS84: GPS Ellipsoid frame
    /// semi-major axis: 6378137.0 [m]
    /// flattening: 1.0/298.2572235630
    WGS84,
    /// WGS72: semi-major axis: 6378135.0 [m]
    /// flattening: 1.0/298.26
    WGS72,
    /// WGS66: semi-major axis: 6378145.0 [m]
    /// flattening: 1.0/298.25
    WGS66,
    /// WGS60: semi-major axis: 6378165.0 [m]
    /// flattening: 1.0/298.3
    WGS60,
    /// PZ90: Glonass Ellipsoid frame
    /// semi-major axis: 6378136.0 [m]
    /// flattening: 1/298.257839303
    PZ90,
    /// BDC, also known as CGCS2000,
    /// is the reference frame used by the
    /// Beidou constellation.
    /// Semi-major axis: 6378137.0 [m]
    /// flattening: 1/298.257222101
    BDC,
    /// GRS80 reference ellipsoid
    /// semi-major axis: 6378137.0 [m]
    /// flattening: 1.0/298.257222101
    GRS80,
    /// Bessel reference ellipsoid
    /// semi-major axis: 6377397.155 [m]
    /// flattening: 1.0/299.1528128
    Bessel,
    /// Airy reference ellipsoid
    /// semi-major axis: 6377563.396 [m]
    /// flattening: 1.0/299.3249646
    Airy,
    /// International reference ellipsoid
    /// semi-major axis: 6378388.0 [m]
    /// flattening: 1.0/297.0
    International,
};

namespace internal {
inline double getMajor(Ellipsoid ellipsoid) {
    switch (ellipsoid) {
        case Ellipsoid::WGS84:
            return 6378137.0;
        case Ellipsoid::WGS72:
            return 6378135.0;
        case Ellipsoid::WGS66:
            return 6378145.0;
        case Ellipsoid::WGS60:
            return 6378165.0;
        case Ellipsoid::PZ90:
            return 6378136.0;
        case Ellipsoid::BDC:
            return 6378137.0;
        case Ellipsoid::GRS80:
            return 6378137.0;
        case Ellipsoid::Bessel:
            return 6377397.155;
        case Ellipsoid::Airy:
            return 6377563.396;
        case Ellipsoid::International:
            return 6378388.0;
        default:
            return 6378137.0;  // unreachable
    }
}

inline double getFlattening(Ellipsoid ellipsoid) {
    switch (ellipsoid) {
        case Ellipsoid::WGS84:
            return 1.0 / 298.257223563;
        case Ellipsoid::WGS72:
            return 1.0 / 298.26;
        case Ellipsoid::WGS66:
            return 1.0 / 298.25;
        case Ellipsoid::WGS60:
            return 1.0 / 298.3;
        case Ellipsoid::PZ90:
            return 1.0 / 298.257839303;
        case Ellipsoid::BDC:
            return 1.0 / 298.257222101;
        case Ellipsoid::GRS80:
            return 1.0 / 298.2572221009;
        case Ellipsoid::Bessel:
            return 299.1528128;
        case Ellipsoid::Airy:
            return 299.3249646;
        case Ellipsoid::International:
            return 297.0;
        default:
            return 1.0 / 298.257223563;  // unreachable
    }
}

inline double getMinor(double major, double flattening) {
    return major * (1.0 - flattening);
}

inline double getSquaredEccentricity(double major, double minor) {
    return ((major * major) - (minor * minor)) / (major * major);
}

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
    auto t = std::cos(lat) * up - std::sin(lat) * nt;
    out_u = std::cos(lon) * t - std::sin(lon) * et;
    out_v = std::sin(lon) * t + std::cos(lon) * et;
    out_w = std::sin(lat) * up + std::cos(lat) * nt;
}

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
    auto t = std::cos(lon) * u + std::sin(lon) * v;
    out_east = -1 * std::sin(lon) * u + std::cos(lon) * v;
    out_north = -1 * std::sin(lat) * t + std::cos(lat) * w;
    out_up = std::cos(lat) * t + std::sin(lat) * w;
}

}  // namespace internal

inline void ecef2geodetic(
    double x,
    double y,
    double z,
    double& out_lat,
    double& out_lon,
    double& out_alt,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double major = internal::getMajor(ellipsoid);
    double flattening = internal::getFlattening(ellipsoid);
    double minor = internal::getMinor(major, flattening);

    double r = std::sqrt(x * x + y * y + z * z);
    double e = std::sqrt(major * major - minor * minor);
    double var = r * r - e * e;
    double u =
        std::sqrt(0.5 * var + 0.5 * std::sqrt(var * var + 4.0 * e * e * z * z));

    double q = std::sqrt(x * x + y * y);
    double hu_e = std::sqrt(u * u + e * e);
    double beta = std::atan(hu_e / u * z / q);

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

inline void aer2enu(
    double az,
    double el,
    double range,
    double& out_e,
    double& out_n,
    double& out_u
) {
    auto r = range * std::cos(el);
    out_e = r * std::sin(az);
    out_n = r * std::cos(az);
    out_u = range * std::sin(el);
}

inline void enu2aer(
    double east,
    double north,
    double up,
    double& out_az,
    double& out_el,
    double& out_range
) {
    double r = std::sqrt(east * east + north * north);
    out_range = std::sqrt(r * r + up * up);
    out_el = std::atan2(up, r);
    const auto pi = (2 * 3.14159265358979311599796346854);

    out_az = std::fmod(pi + std::fmod(std::atan2(east, north), pi), pi);
}

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

inline void geodetic2ecef(
    double lat,
    double lon,
    double alt,
    double& out_x,
    double& out_y,
    double& out_z,
    Ellipsoid ellipsoid = Ellipsoid::WGS84
) {
    double major = internal::getMajor(ellipsoid);
    double flattening = internal::getFlattening(ellipsoid);
    double minor = internal::getMinor(major, flattening);
    double se = internal::getSquaredEccentricity(major, minor);

    double n = major / (std::sqrt(1.0 - se * std::sin(lat) * std::sin(lat)));

    out_x = (n + alt) * std::cos(lat) * std::cos(lon);
    out_y = (n + alt) * std::cos(lat) * std::sin(lon);
    out_z = (n * (minor / major) * (minor / major) + alt) * std::sin(lat);
}

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