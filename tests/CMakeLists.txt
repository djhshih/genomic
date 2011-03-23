
include_directories(${Boost_INCLUDE_DIRS})

set(genomicTest_SRCS genomicTest.cpp ../Sample.cpp ../global.cpp)

add_executable(genomicTest ${genomicTest_SRCS})
target_link_libraries(genomicTest ${Boost_LIBRARIES} -lboost_unit_test_framework)

# Copy data files from source to build tree
file(COPY data/ DESTINATION . FILES_MATCHING REGEX "/.*\\.(cfg|in|ans|cn|seg)$")

add_test(genomicTest genomicTest)