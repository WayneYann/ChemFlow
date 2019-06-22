cmake_minimum_required(VERSION 3.2)

project(counterflow)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories("/usr/include/eigen3")
# include_directories("/home/xzhang/Computation/sundials-4.1.0/instdir/include")
# link_directories("/home/xzhang/Computation/sundials-4.1.0/instdir/lib")

set(SOURCE_FILES main.cpp)

add_executable(counterflow  ${SOURCE_FILES})

# TARGET_LINK_LIBRARIES(counterflow sundials_cvode sundials_nvecserial )
