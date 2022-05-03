# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.12)

project(HESVAT VERSION 1.0 LANGUAGES CXX)

add_compile_options(-O2)
add_compile_options(-std=c++17)
add_compile_options(-march=native)
add_compile_options(-g)

add_executable(test "src/main.cpp" "src/utils.cpp" "src/thread.cpp" "src/testPVAT.cpp" "src/evaluator_Aggregation.cpp"  "src/testSVAT.cpp")

add_executable(testPlain "src/mainPlain.cpp" "src/utils.cpp" "src/testPVAT.cpp")

# Import Microsoft SEAL
find_package(SEAL 3.6.4 EXACT REQUIRED
    PATHS src/external/SEAL-main/cmake 
)
# Import OpenMP 
find_package(OpenMP)
#add_definitions("-fopenmp")

target_link_libraries(test SEAL::seal)
target_link_libraries(test ${OpenMP_CXX_LIBRARIES})

