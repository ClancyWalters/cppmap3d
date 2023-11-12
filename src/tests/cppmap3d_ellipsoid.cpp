#include "doctest/doctest.h"

#include <stdexcept>
#include <vector>
#include "../cppmap3d.hh"
#include "cppmap3d_test_util.hh"

TEST_CASE("flattening") {
    // clang-format off
    std::vector<std::tuple<cppmap3d::Ellipsoid, double>> data_container {
        {cppmap3d::Ellipsoid::Maupertuis, 0.005235602050865236},
        {cppmap3d::Ellipsoid::Plessis, 0.003240020729165458},
        {cppmap3d::Ellipsoid::Everest1830, 0.003324448922118313},
        {cppmap3d::Ellipsoid::Everest1830Modified, 0.003324449295589469},
        {cppmap3d::Ellipsoid::Everest1967, 0.003324449343845343},
        {cppmap3d::Ellipsoid::Airy, 0.00334085067870327},
        {cppmap3d::Ellipsoid::Bessel, 0.0033427731536659813},
        {cppmap3d::Ellipsoid::Clarke1866, 0.0033900753039287908},
        {cppmap3d::Ellipsoid::Clarke1878, 0.003407549790771363},
        {cppmap3d::Ellipsoid::Clarke1860, 0.003407561308111843},
        {cppmap3d::Ellipsoid::Helmert, 0.0033523298109184524},
        {cppmap3d::Ellipsoid::Hayford, 0.003367003387062615},
        {cppmap3d::Ellipsoid::International1924, 0.003367003387062615},
        {cppmap3d::Ellipsoid::Krassovsky1940, 0.0033523298336767685},
        {cppmap3d::Ellipsoid::WGS66, 0.0033528919458556804},
        {cppmap3d::Ellipsoid::Australian, 0.003352891899858333},
        {cppmap3d::Ellipsoid::International1967, 0.003352896192983603},
        {cppmap3d::Ellipsoid::GRS67, 0.0033529237272191623},
        {cppmap3d::Ellipsoid::SA1969, 0.003352891899858333},
        {cppmap3d::Ellipsoid::WGS72, 0.0033527794541680267},
        {cppmap3d::Ellipsoid::GRS80, 0.0033528106811816882},
        {cppmap3d::Ellipsoid::WGS84, 0.0033528106647473664},
        {cppmap3d::Ellipsoid::IERS1989, 0.0033528131102879993},
        {cppmap3d::Ellipsoid::IERS2003, 0.0033528131084554157},
        {cppmap3d::Ellipsoid::Mercury, 0.0009014546199549272},
        {cppmap3d::Ellipsoid::Venus, 0.0},
        {cppmap3d::Ellipsoid::Moon, 0.0012082158679017317},
        {cppmap3d::Ellipsoid::Mars, 0.006123875928193323},
        {cppmap3d::Ellipsoid::Jupyter, 0.06604858798757626},
        {cppmap3d::Ellipsoid::Io, 0.0075968738044488665},
        {cppmap3d::Ellipsoid::Uranus, 0.022927344575296372},
        {cppmap3d::Ellipsoid::Neptune, 0.01708124697141011},
        {cppmap3d::Ellipsoid::Pluto, 0.0},
    };
    // clang-format on

    for (const auto& data : data_container) {
        CAPTURE(data);

        CHECK(
            cppmap3d::internal::getFlattening(std::get<cppmap3d::Ellipsoid>(data
            )) == doctest::Approx(std::get<double>(data))
        );
    }
}

TEST_CASE("ellipsoid") {
    // clang-format off
    std::vector<std::tuple<cppmap3d::Ellipsoid, std::vector<double>>> data_container {
        {cppmap3d::Ellipsoid::Maupertuis,{42.123086280313906, -82.00647850636021, -13462.822154350226}},
        {cppmap3d::Ellipsoid::Plessis,{42.008184833614905, -82.00647850636021, 1566.9219075104988}},
        {cppmap3d::Ellipsoid::Everest1830,{42.01302648557789, -82.00647850636021, 1032.4153744896425}},
        {cppmap3d::Ellipsoid::Everest1830Modified,{42.0130266467127, -82.00647850636021, 1027.7254294115853}},
        {cppmap3d::Ellipsoid::Everest1967,{42.01302648557363, -82.00647850636021, 1033.2243733811288}},
        {cppmap3d::Ellipsoid::Airy,{42.01397060398504, -82.00647850636021, 815.5499438015993}},
        {cppmap3d::Ellipsoid::Bessel,{42.01407537004288, -82.00647850636021, 987.0246149983182}},
        {cppmap3d::Ellipsoid::Clarke1866,{42.01680003414445, -82.00647850636021, 313.90267925120395}},
        {cppmap3d::Ellipsoid::Clarke1878,{42.0177971504227, -82.00647850636021, 380.12002203958457}},
        {cppmap3d::Ellipsoid::Clarke1860,{42.017799612218326, -82.00647850636021, 321.0980872430816}},
        {cppmap3d::Ellipsoid::Helmert,{42.01464497456125, -82.00647850636021, 212.63680219872765}},
        {cppmap3d::Ellipsoid::Hayford,{42.01548834310426, -82.00647850636021, 66.77070154259877}},
        {cppmap3d::Ellipsoid::International1924,{42.01548834310426, -82.00647850636021, 66.77070154259877}},
        {cppmap3d::Ellipsoid::Krassovsky1940,{42.01464632634865, -82.00647850636021, 167.7043859419633}},
        {cppmap3d::Ellipsoid::WGS66,{42.014675415414274, -82.00647850636021, 269.1575142686737}},
        {cppmap3d::Ellipsoid::Australian,{42.01467586302664, -82.00647850636021, 254.17989315657786}},
        {cppmap3d::Ellipsoid::International1967,{42.01467603307557, -82.00647850636021, 256.6883857005818}},
        {cppmap3d::Ellipsoid::GRS67,{42.01467768000789, -82.00647850636021, 254.27066653452297}},
        {cppmap3d::Ellipsoid::SA1969,{42.01467586302664, -82.00647850636021, 254.17989315657786}},
        {cppmap3d::Ellipsoid::WGS72,{42.01466869328149, -82.00647850636021, 278.8216763935984}},
        {cppmap3d::Ellipsoid::GRS80,{42.01467053601299, -82.00647850636021, 276.9137384511387}},
        {cppmap3d::Ellipsoid::WGS84,{42.01467053507479, -82.00647850636021, 276.91369158042767}},
        {cppmap3d::Ellipsoid::WGS84Mean,{41.823366301, -82.0064785, -2.13061272e3}},
        {cppmap3d::Ellipsoid::IERS1989,{42.01467064467172, -82.00647850636021, 277.9191657339711}},
        {cppmap3d::Ellipsoid::IERS2003,{42.01467066257621, -82.00647850636021, 277.320060889772}},
        {cppmap3d::Ellipsoid::Mercury,{41.8430384333997, -82.00647850636021, 3929356.5648451606}},
        {cppmap3d::Ellipsoid::Venus,{41.82336630167669, -82.00647850636021, 317078.15867127385}},
        {cppmap3d::Ellipsoid::Moon,{41.842147614909734, -82.00647850636021, 4631711.995926845}},
        {cppmap3d::Ellipsoid::Mars,{42.00945156056578, -82.00647850636021, 2981246.073616111}},
        {cppmap3d::Ellipsoid::Jupyter,{75.3013267078341, -82.00647850636021, -61782040.202975556}},
        {cppmap3d::Ellipsoid::Io,{41.82422244977044, -82.00647850636021, 6367054.626528843}},
        {cppmap3d::Ellipsoid::Uranus,{47.69837228395133, -82.00647850636021, -18904824.4361074}},
        {cppmap3d::Ellipsoid::Neptune,{45.931317431546425, -82.00647850636021, -18194050.781948525}},
        {cppmap3d::Ellipsoid::Pluto,{41.82336630167669, -82.00647850636021, 5180878.158671274}}
    };
    // clang-format on

    for (const auto& data : data_container) {
        CAPTURE(data);
        double x = 660e3, y = -4700e3, z = 4247e3;
        double lat, lon, alt;
        cppmap3d::ecef2geodetic(
            x,
            y,
            z,
            lat,
            lon,
            alt,
            std::get<cppmap3d::Ellipsoid>(data)
        );

        CHECK(degrees(lat) == doctest::Approx(std::get<1>(data)[0]));
        CHECK(degrees(lon) == doctest::Approx(std::get<1>(data)[1]));
        CHECK(alt == doctest::Approx(std::get<1>(data)[2]));
    }
}