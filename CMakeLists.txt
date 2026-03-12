cmake_minimum_required(VERSION 3.20)

project(math_engine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()    


# Library target
file(GLOB_RECURSE MATH_ENGINE_SOURCES CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/src/*.cpp
)

add_library(math_engine ${MATH_ENGINE_SOURCES})

target_include_directories(math_engine
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)


# Testing
enable_testing()
add_subdirectory(tests)