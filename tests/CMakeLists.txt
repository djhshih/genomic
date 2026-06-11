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

file(COPY data/ DESTINATION . FILES_MATCHING REGEX "/.*\\.(cfg|in|ans|cn|seg|cnas|segas|txt|tsv)$")
file(COPY cbs_compare.R DESTINATION .)
file(COPY nclist_compare.R DESTINATION .)
file(COPY data/cbs_case1_perm_alt_expected.tsv DESTINATION .)
file(COPY data/cbs_case1_hybrid_expected.tsv DESTINATION .)
file(COPY data/cbs_case1_hybrid_alt_expected.tsv DESTINATION .)
file(COPY data/cbs_case2_weighted_perm_alt_expected.tsv DESTINATION .)
file(COPY data/cbs_case2_weighted_hybrid_expected.tsv DESTINATION .)
file(COPY data/cbs_case2_weighted_hybrid_alt_expected.tsv DESTINATION .)
file(COPY data/cbs_case3_noisy_perm_alt_expected.tsv DESTINATION .)
file(COPY data/cbs_case3_noisy_hybrid_expected.tsv DESTINATION .)
file(COPY data/cbs_case3_noisy_hybrid_alt_expected.tsv DESTINATION .)
file(COPY data/cbs_case4_noisy_perm_alt_expected.tsv DESTINATION .)
file(COPY data/cbs_case4_noisy_hybrid_expected.tsv DESTINATION .)
file(COPY data/cbs_case4_noisy_hybrid_alt_expected.tsv DESTINATION .)

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
