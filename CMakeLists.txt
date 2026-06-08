cmake_minimum_required(VERSION 4.0)
project(cna)
enable_testing()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(CNA_ENABLE_WARNINGS "Enable extra compiler warnings" ON)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	message("Building Debug version")
	set(debug_build 1)
else()
	set(debug_build 0)
endif()

set(version_major 0)
set(version_minor 0)

configure_file(
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_BINARY_DIR}/config.h"
)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules" ${CMAKE_MODULE_PATH})

find_package(Boost REQUIRED CONFIG)

set(lib_sources
	src/cna_common.cpp
	src/cna_clean.cpp
	src/cna_convert.cpp
	src/cna_filter.cpp
	lib/global.cpp
	lib/logging.cpp
	lib/parse.cpp
	lib/Sample.cpp
	lib/SampleSet.cpp
	lib/GenericSampleSet.cpp
	lib/Marker.cpp
)

add_library(cna_lib ${lib_sources})
target_include_directories(cna_lib PUBLIC ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src)
set_target_properties(cna_lib PROPERTIES OUTPUT_NAME "cna")

if(CNA_ENABLE_WARNINGS)
	target_compile_options(cna_lib PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
	)
endif()

add_executable(cna src/cna.cpp)
target_include_directories(cna PRIVATE ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(cna cna_lib boost_program_options)

if(CNA_ENABLE_WARNINGS)
	target_compile_options(cna PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
	)
endif()

install(TARGETS cna RUNTIME DESTINATION bin)

add_subdirectory(tests)
#add_subdirectory(examples)
