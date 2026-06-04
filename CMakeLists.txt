cmake_minimum_required(VERSION 3.5)
project(genomic)
enable_testing()

# configure project properties

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

# include custom modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules" ${CMAKE_MODULE_PATH})



find_package(Boost REQUIRED)

include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${BOOST_INCLUDES})

set(genomic_LIB_SRCS
	src/genomic_common.cpp
	lib/global.cpp
	lib/Sample.cpp
	lib/SampleSet.cpp
	lib/GenericSampleSet.cpp
	lib/Marker.cpp
)

add_library(genomic_lib ${genomic_LIB_SRCS})

add_executable(genomic src/genomic.cpp)
target_link_libraries(genomic genomic_lib boost_program_options)

install(TARGETS genomic RUNTIME DESTINATION bin)

add_subdirectory(tests)
#add_subdirectory(examples)
