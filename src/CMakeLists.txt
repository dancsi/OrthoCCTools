include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(external)
add_subdirectory(common)
add_subdirectory(scoring)
add_subdirectory(fastscore)
add_subdirectory(pyccscore)
add_subdirectory(solver)
add_subdirectory(setbuilder)