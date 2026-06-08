cmake_minimum_required(VERSION 4.0)

set(CMAKE_BUILD_TYPE Debug)

find_package(Boost REQUIRED CONFIG)

add_executable(genomic_test genomic_test.cpp)
target_include_directories(genomic_test PRIVATE ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
target_link_libraries(genomic_test genomic_lib boost_unit_test_framework)

if(GENOMIC_ENABLE_WARNINGS)
	target_compile_options(genomic_test PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
	)
endif()

# Copy data files from source to build tree
file(COPY data/ DESTINATION . FILES_MATCHING REGEX "/.*\\.(cfg|in|ans|cn|seg|cnas|segas|txt)$")

add_test(genomic_test genomic_test)

add_executable(genomic_cn_test genomic_cn_test.cpp)
target_include_directories(genomic_cn_test PRIVATE ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
target_link_libraries(genomic_cn_test genomic_lib boost_unit_test_framework)

if(GENOMIC_ENABLE_WARNINGS)
	target_compile_options(genomic_cn_test PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
	)
endif()

add_test(genomic_cn_test genomic_cn_test)
