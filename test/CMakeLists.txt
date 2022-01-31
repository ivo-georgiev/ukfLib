cmake_minimum_required(VERSION 3.3)
project(kalman)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
		"${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

	add_compile_options(-Wall -Wextra -Werror)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

	add_compile_options(/W4 /WX /EHsc)
endif()

add_compile_options(-std=c99 --pedantic)
include_directories(../cfg ../lib)
add_executable (kalman ../cfg/ukfCfg.c ../cfg/ukfCfg1.c ../lib/ukfLib.c ukfTest.c ../lib/mtxLib.c)
target_link_libraries(kalman m)


enable_testing()
add_test( testname kalman )

