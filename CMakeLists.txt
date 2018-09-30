cmake_minimum_required( VERSION 3.2 )
project( SCPT LANGUAGES CXX )

set(CMAKE_CXX_STANDARD 17)

add_library( stlab INTERFACE )
target_include_directories( stlab INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/dependecies/libraries/stlab/concurrency/ )

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/src/SCPT/SCol/ )
add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/src/SCPT/Threaded/ )

