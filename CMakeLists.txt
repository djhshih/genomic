project(genomic)
cmake_minimum_required(VERSION 2.8)
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

include_directories(${PROJECT_BINARY_DIR} ${BOOST_INCLUDES})

set(genomic_SRCS genomic.cpp genomic_common.cpp global.cpp SampleSet.cpp GenericSampleSet.cpp Marker.cpp)

add_executable(genomic ${genomic_SRCS})
target_link_libraries(genomic boost_program_options)

install(TARGETS genomic RUNTIME DESTINATION bin)

add_subdirectory(tests)
add_subdirectory(examples)
add_subdirectory(graphic)
add_subdirectory(command)
