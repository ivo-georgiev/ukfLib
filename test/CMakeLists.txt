cmake_minimum_required(VERSION 3.3)
project(kalman)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
		"${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

	set(warnings "-Wall -Wextra -Werror")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

	set(warnings "/W4 /WX /EHsc")
endif()
if (NOT CONFIGURED_ONCE)
	set(CMAKE_CXX_FLAGS "${warnings}"
		CACHE STRING "Flags used by the compiler during all build types." FORCE)

	set(CMAKE_C_FLAGS
		"${warnings}"
		CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c99")
add_executable (kalman ../cfg/ukfCfg.c ../cfg/ukfCfg1.c ../lib/ukfLib.c ukfTest.c ../lib/mtxLib.c)
target_link_libraries(kalman m)


enable_testing()
add_test( testname kalman )

