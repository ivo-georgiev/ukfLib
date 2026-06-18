cmake_minimum_required(VERSION 3.3)
project(kalman)

# Compiler configuration for strict warnings
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Werror)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    add_compile_options(/W4 /WX /EHsc)
endif()

# Set standard C99 build flags
add_compile_options(-std=c99 --pedantic)

# Add the new submodule directory to include paths so mtxLib.h can be found
include_directories(.. ../cfg ../lib ../lib/mtxLib)


# Build the executable with the updated submodule path for mtxLib.c
add_executable (kalman 
    ../cfg/ukfCfg.c 
    ../cfg/ukfCfg1.c 
    ../lib/ukfLib.c 
    ukfTest.c 
    ../lib/mtxLib/mtxLib.c
)

# Link the standard math library
target_link_libraries(kalman m)

# Enable testing framework
enable_testing()
add_test(testname kalman)
