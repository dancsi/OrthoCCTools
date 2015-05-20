cmake_minimum_required(VERSION 3.1)
project(source_MaxCliquePara_v2_2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -I.")

set(SOURCE_FILES
    BB_ColorRSort.h
    BB_GreedyColorSort.h
    BitSet.cpp
    BitSet.h
    DegreeAndNumberSort.h
    DegreeSort.h
    Dimacs.cpp
    Dimacs.h
    GreedyColorSort.h
    KillTimer.h
    main.cpp
    MaximumCliqueBase.h
    McrBB.h
    ParallelMaximumClique.h
    SteadyGreedyColorSort.h
    Timer.h)

add_executable(maxClique22 ${SOURCE_FILES})