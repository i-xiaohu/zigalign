cmake_minimum_required(VERSION 3.17)
project(zigalign)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")

add_executable(zigalign
        kseq.h
        input_fasta.h
        input_fasta.cpp
        zigalign.h
        zigalign_str.cpp
        zigalign_etr.cpp
        main.cpp)
find_package(OpenMP REQUIRED)
target_link_libraries(zigalign -lz OpenMP::OpenMP_CXX)