set(FASTSCORE_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/fastscore.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEngine.cpp
)

set(FASTSCORE_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEngine.h
)

set(FASTSCORE_ALL ${FASTSCORE_SOURCES} ${FASTSCORE_HEADERS})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scores.dat ${CMAKE_CURRENT_BINARY_DIR}/scores.dat COPYONLY)

add_executable(fastscore ${FASTSCORE_ALL} ${COMMON_ALL})
target_link_libraries(fastscore spdlog flags MemoryMapped ${EXTERNAL_LIBRARIES})