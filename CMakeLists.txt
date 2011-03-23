project(genomic)
cmake_minimum_required(VERSION 2.8)
enable_testing()

find_package(Boost REQUIRED)

include_directories(${BOOST_INCLUDES})

set(genomic_SRCS genomic.cpp Sample.cpp global.cpp)

add_executable(genomic ${genomic_SRCS})
add_subdirectory(tests)