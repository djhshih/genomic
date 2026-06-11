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

execute_process(
	COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/nclist_generate.R
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	RESULT_VARIABLE NCLIST_COMPARE_RESULT
	OUTPUT_QUIET
	ERROR_QUIET
)
if(NOT NCLIST_COMPARE_RESULT EQUAL 0)
	message(WARNING "Failed to run tests/nclist_compare.R")
endif()

execute_process(
	COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/cbs_generate.R
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	RESULT_VARIABLE CBS_COMPARE_RESULT
	OUTPUT_QUIET
	ERROR_QUIET
)
if(NOT CBS_COMPARE_RESULT EQUAL 0)
	message(WARNING "Failed to run tests/cbs_compare.R")
endif()

execute_process(
	COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/smooth_generate.R
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	RESULT_VARIABLE SMOOTH_CNA_GENERATE_RESULT
	OUTPUT_QUIET
	ERROR_QUIET
)
if(NOT SMOOTH_CNA_GENERATE_RESULT EQUAL 0)
	message(WARNING "Failed to run tests/smooth_generate.R")
endif()

file(COPY data/ DESTINATION . FILES_MATCHING REGEX "/.*\\.(cfg|in|ans|cn|seg|cnas|segas|txt|tsv)$")

add_test(cna_test cna_test)

add_executable(cbs_test cbs_test.cpp)
target_include_directories(cbs_test PRIVATE ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
target_link_libraries(cbs_test cna_lib boost_unit_test_framework)

if(CNA_ENABLE_WARNINGS)
	target_compile_options(cbs_test PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
	)
endif()

add_test(cbs_test cbs_test)

add_executable(smooth_test smooth_test.cpp)
target_include_directories(smooth_test PRIVATE ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
target_link_libraries(smooth_test cna_lib boost_unit_test_framework)

if(CNA_ENABLE_WARNINGS)
	target_compile_options(smooth_test PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
	)
endif()

add_test(smooth_test smooth_test)

add_executable(nclist_test nclist_test.cpp)
target_include_directories(nclist_test PRIVATE ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
target_link_libraries(nclist_test cna_lib boost_unit_test_framework)

if(CNA_ENABLE_WARNINGS)
	target_compile_options(nclist_test PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
	)
endif()

add_test(nclist_test nclist_test)

#add_executable(cna_cn_test cna_cn_test.cpp)
# target_include_directories(cna_cn_test PRIVATE ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
# target_link_libraries(cna_cn_test cna_lib boost_unit_test_framework)

# if(CNA_ENABLE_WARNINGS)
# 	target_compile_options(cna_cn_test PRIVATE
# 		$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
# 	)
# endif()

#add_test(cna_cn_test cna_cn_test)
