cmake_minimum_required(VERSION 2.8.4)
project(test)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
set(GTEST_ROOT ${CMAKE_SOURCE_DIR}/gtest-1.7.0)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

include_directories(
    ${GTEST_ROOT}/include
    ${CMAKE_SOURCE_DIR}/include
)

add_executable(gtest_pixel
    PixelTest.cpp
)
target_link_libraries(gtest_pixel
    ${CMAKE_SOURCE_DIR}/lib/libimage.a
    ${GTEST_ROOT}/lib/libgtest.a
    ${GTEST_ROOT}/lib/libgtest_main.a
)