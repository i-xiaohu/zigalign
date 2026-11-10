cmake_minimum_required(VERSION 3.17)
project(zigalign)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3")

add_executable(zigalign
        edlib.cpp
        edlib.h
        kseq.h
        ksw2.h
        ksw2_gg2.c
        input_fasta.h
        input_fasta.cpp
        zigalign.h
        zigalign_str.cpp
        zigalign_etr.cpp
        main.cpp
        wraparound.cpp
        wraparound.h
        union_set.h)
find_package(OpenMP REQUIRED)
target_link_libraries(zigalign -lz OpenMP::OpenMP_CXX)