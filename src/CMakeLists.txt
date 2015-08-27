cmake_minimum_required(VERSION 3.2)
project(CCOrtho)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    mcqd_para/BB_ColorRSort.h
    mcqd_para/BB_GreedyColorSort.h
    mcqd_para/BitSet.cpp
    mcqd_para/BitSet.h
    mcqd_para/DegreeAndNumberSort.h
    mcqd_para/DegreeSort.h
    mcqd_para/GreedyColorSort.h
    mcqd_para/KillTimer.h
    mcqd_para/MaximumCliqueBase.h
    mcqd_para/McrBB.h
    mcqd_para/ParallelMaximumClique.h
    mcqd_para/SteadyGreedyColorSort.h
    mcqd_para/Timer.h
    options.h
    solver.cpp)

set(FASTSCORE_SOURCES
    fastscore/Interaction.cpp
    fastscore/Interaction.h
    fastscore/main.cpp)

add_executable(solver ${SOURCE_FILES})
add_executable(fastscore ${FASTSCORE_SOURCES})
