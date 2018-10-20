cmake_minimum_required( VERSION 3.2 )
project( SCPT LANGUAGES CXX )

set(CMAKE_CXX_STANDARD 17)

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/range-v3 )

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/src/SCPT/SCol/ )
add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/src/SCPT/Simple_Collision/ )


