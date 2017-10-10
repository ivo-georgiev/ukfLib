cmake_minimum_required(VERSION 3.3)
project(kalman)

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c99")

add_executable (kalman ukfCfg.c ukfLib.c test2.c mtxLib.c)
target_link_libraries(kalman m)


enable_testing()
add_test( testname kalman )

