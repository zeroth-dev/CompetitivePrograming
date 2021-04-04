cmake_minimum_required( VERSION 3.16)

project(Rosalind)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lcurlpp -lcurl")
add_compile_options("-lcurlpp -lcurl")
add_executable(Rosalind Source.cpp)
TARGET_LINK_LIBRARIES(Rosalind -lcurl -lcurlpp )
