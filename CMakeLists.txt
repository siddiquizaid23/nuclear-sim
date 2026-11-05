cmake_minimum_required(VERSION 3.16)
project(nuclear_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(nuclear_sim src/main.cpp 
    src/reactor_physics.cpp)
    target_compile_options(nuclear_sim PRIVATE -Wall -Wextra)
target_include_directories(nuclear_sim PRIVATE include)
