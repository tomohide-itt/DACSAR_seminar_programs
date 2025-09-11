cmake_minimum_required(VERSION 3.18)

project( cmake_ex CXX )

add_subdirectory( sum )

add_executable( main main.cpp )

target_link_libraries( main PRIVATE sum )


