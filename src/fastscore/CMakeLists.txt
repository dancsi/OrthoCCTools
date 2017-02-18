set(FASTSCORE_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/fastscore.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEnginePotapov.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEngineBCIPA.cpp
)

set(FASTSCORE_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringHelper.h
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEnginePotapov.h
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEngineBCIPA.h
)

set(FASTSCORE_ALL ${FASTSCORE_SOURCES} ${FASTSCORE_HEADERS})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scores.dat ${CMAKE_CURRENT_BINARY_DIR}/scores.dat COPYONLY)

add_executable(fastscore ${FASTSCORE_ALL} ${COMMON_ALL})
target_link_libraries(fastscore spdlog flags MemoryMapped ${EXTERNAL_LIBRARIES})