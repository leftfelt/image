cmake_minimum_required(VERSION 2.8.4)
project(test)

set(GTEST_ROOT ${CMAKE_SOURCE_DIR}/gtest-1.7.0)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_definitions("-std=c++11")

include_directories(
    ${GTEST_ROOT}/include
    ${CMAKE_SOURCE_DIR}/include
)

add_executable(gtest_image
    PixelTest.cpp
    ImageTest.cpp
)

target_link_libraries(gtest_image
    ${CMAKE_SOURCE_DIR}/lib/libimage.a
    ${GTEST_ROOT}/lib/libgtest.a
    ${GTEST_ROOT}/lib/libgtest_main.a
)