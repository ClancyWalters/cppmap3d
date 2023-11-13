.. cppmap3d documentation master file, created by
   sphinx-quickstart on Sun Nov 12 22:48:30 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

cppmap3d Documentation
====================================

.. toctree::
   :maxdepth: 2

https://github.com/ClancyWalters/cppmap3d

Implementation details and examples
=========================================
This library was written as a baseline for 3D coordinate conversions. It's primarily
a direct adaptation of pymap3d and so will share performance characteristics of that
libraries' choice of algorithms. It's intended to have high compatability, even with old 
compilers, and avoids sacrificing performance where possible. 

The library is intended to be minimally invasive pulling in only <cmath> from the STL and
using no data structures. Compatability and avoiding the STL result in the heavy use of 
out paramters, if this is undesirable, I recommend writing a wrapper using structs or std::tuple<>.

.. code-block:: cpp

      #import "cppmap3d.hh"
      
      double x, y, z;
      double lat, lon, alt;
      double az, el, range;
      cppmap3d::ecef2aer(x, y, z, lat, lon, alt, az, el, range);


One concession to performance is a #define exists for using a significantly faster
ecef2geodetic. This breaks parity with pymap3d.

.. code-block:: cpp

      #define CPPMAP3D_ECEF2GEODETIC_OLSON
      #import "cppmap3d.hh"
      
      double x, y, z;
      double lat, lon, alt;
      cppmap3d::ecef2geodetic(x, y, z, lat, lon, alt); //roughly 2x speed of You (2000)



Development / Running DocTest & Nanobench
=========================================
Clone the project recursively

.. code-block:: console

      git clone --recursive https://github.com/ClancyWalters/cppmap3d.git


Setup vcpkg

.. code-block:: console

      .\vcpkg\bootstrap-vcpkg.bat
      .\vcpkg\vcpkg.exe install
      .\vcpkg\vcpkg.exe integrate install
      

Build the project, tests, and Docs. Obviously choose an appropriate present.

.. code-block:: console

      cmake --build --preset=windows-x64-debug

Run the tests

.. code-block:: console

      .\out\build\windows-x64-debug\cppmap3d_tests.exe

Run the benchmarking (should be done in release)

.. code-block:: console

      .\out\build\windows-x64-release\cppmap3d_benchmarks.exe

documentation is compiled locally at out/build/windows-x64-release/docs/sphinx/index.html

Docs
====

.. doxygenfile:: cppmap3d.hh