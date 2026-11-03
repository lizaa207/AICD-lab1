cmake_minimum_required(VERSION 3.10)
project(set_lab CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(set_lab main.cpp)

if(MSVC)
  target_compile_options(set_lab PRIVATE /W4)
else()
  target_compile_options(set_lab PRIVATE -Wall -Wconversion -Wextra -Wpedantic)
endif()