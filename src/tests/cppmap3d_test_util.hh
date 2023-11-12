#ifndef CPPMAP3D_TEST_UTIL_GUARD
#define CPPMAP3D_TEST_UTIL_GUARD

#include "doctest/doctest.h"

#include <corecrt_math_defines.h>
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <vector>

inline double radians(double degrees) {
    return degrees * (M_PI / 180.0);
}

inline double degrees(double radians) {
    return radians * (180.0 / M_PI);
}

// Providing a hint to doctest to allow proper printing of input values
namespace doctest {
template <>
struct StringMaker<std::vector<std::vector<double>>> {
    static String convert(const std::vector<std::vector<double>>& arguments) {
        doctest::String out = "{";
        for (const auto& argument_group : arguments) {
            out += "{";
            for (const auto& argument : argument_group) {
                out += doctest::toString(argument);
                if (&argument != &argument_group.back()) {
                    out += ", ";
                }
            }
            out += "}";
            if (&argument_group != &arguments.back()) {
                out += ", ";
            }
        }
        out += "}";
        return out;
    }
};
}  // namespace doctest

#endif