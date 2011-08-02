project(genomic)
cmake_minimum_required(VERSION 2.8)
enable_testing()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	message("Building Debug version")
	set(genomic_DEBUG 1)
else()
	set(genomic_DEBUG 0)
endif()

set(genomic_VERSION_MAJOR 0)
set(genomic_VERSION_MINOR 0)

configure_file(
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_BINARY_DIR}/config.h"
)


find_package(Boost REQUIRED)

include_directories(${PROJECT_BINARY_DIR} ${BOOST_INCLUDES})

set(genomic_SRCS genomic.cpp global.cpp Sample.cpp)

add_executable(genomic ${genomic_SRCS})
target_link_libraries(genomic boost_program_options)

add_subdirectory(tests)
add_subdirectory(examples)