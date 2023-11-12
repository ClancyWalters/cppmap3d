#ifndef CPPMAP3D_TEST_UTIL_GUARD
#define CPPMAP3D_TEST_UTIL_GUARD

#include "doctest/doctest.h"

#include <corecrt_math_defines.h>
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <tuple>
#include <vector>

#include "../cppmap3d.hh"

inline double radians(double degrees) {
    return degrees * (M_PI / 180.0);
}

inline double degrees(double radians) {
    return radians * (180.0 / M_PI);
}

// Providing a hint to doctest to allow proper printing of input values
namespace doctest {
template <typename T>
struct StringMaker<std::vector<T>> {
    static String convert(const std::vector<T>& arguments) {
        doctest::String out = "{";

        for (const T& argument : arguments) {
            out += doctest::toString(argument);
            if (&argument != &arguments.back()) {
                out += ", ";
            }
        }

        out += "}";
        return out;
    }
};

template <>
struct StringMaker<std::vector<double>> {
    static String convert(const std::vector<double>& arguments) {
        doctest::String out = "{";

        for (const auto& argument : arguments) {
            out += doctest::toString(argument);
            if (&argument != &arguments.back()) {
                out += ", ";
            }
        }

        out += "}";
        return out;
    }
};

template <>
struct StringMaker<std::tuple<cppmap3d::Ellipsoid, double>> {
    static String convert(
        const std::tuple<cppmap3d::Ellipsoid, double>& argument
    ) {
        doctest::String out = "{";

        out += doctest::toString(std::get<cppmap3d::Ellipsoid>(argument));
        out += ", ";
        out += doctest::toString(std::get<double>(argument));

        out += "}";
        return out;
    }
};

template <>
struct StringMaker<std::tuple<cppmap3d::Ellipsoid, std::vector<double>>> {
    static String convert(
        const std::tuple<cppmap3d::Ellipsoid, std::vector<double>>& argument
    ) {
        doctest::String out = "{";

        out += doctest::toString(std::get<cppmap3d::Ellipsoid>(argument));
        out += ", ";
        out += doctest::toString(std::get<std::vector<double>>(argument));

        out += "}";
        return out;
    }
};

}  // namespace doctest

#endif