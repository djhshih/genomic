
include_directories(${Boost_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

set(genomicTest_SRCS genomicTest.cpp ../global.cpp ../Sample.cpp)

add_executable(genomicTest ${genomicTest_SRCS})
target_link_libraries(genomicTest boost_unit_test_framework)

# Copy data files from source to build tree
file(COPY data/ DESTINATION . FILES_MATCHING REGEX "/.*\\.(cfg|in|ans|cn|seg)$")

add_test(genomicTest genomicTest)