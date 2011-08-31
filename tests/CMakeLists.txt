set(CMAKE_BUILD_TYPE Debug)

find_package(Boost REQUIRED)

include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
include_directories(${Boost_INCLUDE_DIRS})

set(genomic_test_SRCS genomic_test.cpp ../global.cpp ../SampleSet.cpp ../GenericSampleSet.cpp ../Marker.cpp)

add_executable(genomic_test ${genomic_test_SRCS})
target_link_libraries(genomic_test boost_unit_test_framework)

# Copy data files from source to build tree
file(COPY data/ DESTINATION . FILES_MATCHING REGEX "/.*\\.(cfg|in|ans|cn|seg|cnas|segas)$")

add_test(genomic_test genomic_test)
