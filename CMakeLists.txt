cmake_minimum_required(VERSION 4.0)
project(genomic)
enable_testing()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(GENOMIC_ENABLE_WARNINGS "Enable extra compiler warnings" ON)

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



find_package(Boost REQUIRED CONFIG)

set(genomic_LIB_SRCS
	src/genomic_common.cpp
	lib/global.cpp
	lib/logging.cpp
	lib/Sample.cpp
	lib/SampleSet.cpp
	lib/GenericSampleSet.cpp
	lib/Marker.cpp
)

add_library(genomic_lib ${genomic_LIB_SRCS})
target_include_directories(genomic_lib PUBLIC ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src)
set_target_properties(genomic_lib PROPERTIES OUTPUT_NAME "genomic")

if(GENOMIC_ENABLE_WARNINGS)
	target_compile_options(genomic_lib PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
	)
endif()

add_executable(genomic src/genomic.cpp)
target_include_directories(genomic PRIVATE ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(genomic genomic_lib boost_program_options)

if(GENOMIC_ENABLE_WARNINGS)
	target_compile_options(genomic PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
	)
endif()

install(TARGETS genomic RUNTIME DESTINATION bin)


add_subdirectory(tests)
#add_subdirectory(examples)
