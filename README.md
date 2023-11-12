# cppmap3d
Header only single file C++ Implementation of pymap3d. API is similar to the $1000 Matlab Mapping toolbox.

# Similar toolboxes in other languages

* [Matlab, GNU Octave](https://github.com/geospace-code/matmap3d)
* [Fortran](https://github.com/geospace-code/maptran3d)
* [Rust](https://github.com/gberrante/map_3d)
* [Python](https://github.com/geospace-code/pymap3d/tree/main)

## Usage
```C++
#include "cppmap3d.hh"

double x, y, z;
cppmap3d::geodetic2ecef(lat, lon, alt, x, y, z);

double az, el, range;
cppmap3d::geodetic2aer(lat, lon, alt, observer_lat, observer_lon, observer_alt, az, el, range);
```

If speed is important the library includes an alternative ecef2geodetic implementation that is faster according to [this comparison](https://github.com/planet36/ecef-geodetic/tree/main). A microbenchmark is included in `src/benchmarks` indicating about 2x the performance. However, this implementation does not handle degenerate cases as well.

```C++
#define CPPMAP3D_ECEF2GEODETIC_OLSON
#include "cppmap3d.hh"

double lat, lon, alt;
cppmap3d::ecef2geodetic(x, y, z, lat, lon, alt);
```

## Install

Copy `cppmap3d.hh` into your project files.

## Prerequisites

Tests run at C++11, use at own risk for older compiler versions.

## Functions

```
aer2ecef  aer2enu  aer2geodetic  aer2ned
ecef2aer  ecef2enu  ecef2geodetic  ecef2ned
enu2aer  enu2ecef   enu2geodetic
geodetic2aer  geodetic2ecef  geodetic2enu  geodetic2ned
ned2aer  ned2ecef   ned2geodetic
```
### Abbreviations:

* [AER: Azimuth, Elevation, Range](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
* [ECEF: Earth-centered, Earth-fixed](https://en.wikipedia.org/wiki/Earth-centered,_Earth-fixed_coordinate_system)
* [ENU: East North Up](https://en.wikipedia.org/wiki/Axes_conventions#Ground_reference_frames:_ENU_and_NED)
* [NED: North East Down](https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates)

## Limitations

Unlike pymap3d this library cannot take advantage of SIMD by processing large vectorized inputs. This is non-trivial in C++ but feel free to submit a PR.

Doubles are used for all operations, if anyone wants to refactor for arbitrary precision, submit a PR.