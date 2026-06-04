cmake_minimum_required(VERSION 3.5)

set(CMAKE_BUILD_TYPE Debug)

find_package(Boost REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
include_directories(${Boost_INCLUDE_DIRS})

add_executable(genomic_test genomic_test.cpp)
target_link_libraries(genomic_test genomic_lib boost_unit_test_framework)

# Copy data files from source to build tree
file(COPY data/ DESTINATION . FILES_MATCHING REGEX "/.*\\.(cfg|in|ans|cn|seg|cnas|segas)$")

add_test(genomic_test genomic_test)
