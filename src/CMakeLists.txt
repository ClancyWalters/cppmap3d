# FetchContent protection
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL "${CMAKE_SOURCE_DIR}/src")
	add_subdirectory( benchmarks )
	add_subdirectory( tests )
endif()

# Create a cppmap3d target for FetchContent
add_library(cppmap3d INTERFACE)
add_library(cppmap3d::cppmap3d ALIAS cppmap3d)

target_include_directories(cppmap3d INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)