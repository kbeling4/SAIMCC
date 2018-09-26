add_executable( Scol Scol.cpp)

find_package(Boost 1.45.0) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    target_link_libraries(Scol ${Boost_LIBRARIES})
endif()
