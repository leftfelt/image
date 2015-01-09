cmake_minimum_required(VERSION 2.8.4)
project(image)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/library)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/library)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/library)
set(image_source
    source/Area.cpp
    source/Pixel.cpp
    source/Point.cpp
    source/Math.cpp
    source/Image.cpp
    source/ImageIO.cpp
    source/ImageUtil.cpp
)
include_directories(include)

add_library(image SHARED ${image_source})
add_library(image-static STATIC ${image_source})
set_target_properties(image-static PROPERTIES OUTPUT_NAME "image")

add_library(image-win SHARED ${image_source})
set_target_properties(image-win PROPERTIES OUTPUT_NAME "image")
set_target_properties(image-win PROPERTIES PREFIX "")
set_target_properties(image-win PROPERTIES SUFFIX ".dll")

add_library(image-win-static STATIC ${image_source})
set_target_properties(image-win-static PROPERTIES OUTPUT_NAME "image")
set_target_properties(image-win-static PROPERTIES PREFIX "")
set_target_properties(image-win-static PROPERTIES SUFFIX ".lib")