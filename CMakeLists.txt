cmake_minimum_required(VERSION 3.20)
project(Base256 VERSION 4.6.5 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(base256 INTERFACE)
add_library(Base256::base256 ALIAS base256)
target_include_directories(base256 INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

include(CTest)
if(BUILD_TESTING)
    add_executable(test_float256 tests/test_float256.cpp)
    target_link_libraries(test_float256 PRIVATE base256)
    add_executable(test_extra_math tests/test_extra_math.cpp)
    target_link_libraries(test_extra_math PRIVATE base256)
    if(MSVC)
        target_compile_options(test_float256 PRIVATE /W4 /permissive-)
        target_compile_options(test_extra_math PRIVATE /W4 /permissive-)
    else()
        target_compile_options(test_float256 PRIVATE -Wall -Wextra -Wpedantic)
        target_compile_options(test_extra_math PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    add_test(NAME float256 COMMAND test_float256)
    add_test(NAME extra_math COMMAND test_extra_math)
endif()

add_executable(example_basic examples/basic.cpp)
target_link_libraries(example_basic PRIVATE base256)

add_executable(example_high_precision examples/high_precision.cpp)
target_link_libraries(example_high_precision PRIVATE base256)

install(DIRECTORY include/base256 DESTINATION include)
