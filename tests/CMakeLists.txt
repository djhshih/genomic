cmake_minimum_required(VERSION 4.0)

set(CMAKE_BUILD_TYPE Debug)

find_package(Boost REQUIRED CONFIG)

add_executable(cna_test cna_test.cpp)
target_include_directories(cna_test PRIVATE ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
target_link_libraries(cna_test cna_lib boost_unit_test_framework)

if(CNA_ENABLE_WARNINGS)
	target_compile_options(cna_test PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
	)
endif()

file(COPY data/ DESTINATION . FILES_MATCHING REGEX "/.*\\.(cfg|in|ans|cn|seg|cnas|segas|txt)$")

add_test(cna_test cna_test)

#add_executable(cna_cn_test cna_cn_test.cpp)
# target_include_directories(cna_cn_test PRIVATE ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
# target_link_libraries(cna_cn_test cna_lib boost_unit_test_framework)

# if(CNA_ENABLE_WARNINGS)
# 	target_compile_options(cna_cn_test PRIVATE
# 		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
# 	)
# endif()

#add_test(cna_cn_test cna_cn_test)
