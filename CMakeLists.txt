cmake_minimum_required( VERSION 3.26.4 )

# Export compile commands for clangd
set( CMAKE_EXPORT_COMPILE_COMMANDS TRUE )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

option( SANITIZE "Enable selected sanitizer in build" OFF )

project( "cppmap3d" )

add_subdirectory( src )

# Copy compile commands for clangd
add_custom_target(
	copy-compile-commands ALL
	${CMAKE_COMMAND} -E copy_if_different
	${CMAKE_BINARY_DIR}/compile_commands.json
	${CMAKE_CURRENT_LIST_DIR}
)