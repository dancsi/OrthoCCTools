set(SCORING_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEnginePotapov.cpp
)

set(SCORING_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/CIPAHelper.h
	${CMAKE_CURRENT_SOURCE_DIR}/ICIPAHelper.h
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEngineBCIPA.h
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEngineQCIPA.h
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEngineICIPA.h
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringEnginePotapov.h
	${CMAKE_CURRENT_SOURCE_DIR}/ScoringHelper.h
)

set(SCORING_ALL ${SCORING_SOURCES} ${SCORING_HEADERS})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scores.dat ${CMAKE_CURRENT_BINARY_DIR}/../scores.dat COPYONLY)

add_library(scoring ${SCORING_ALL})

target_include_directories(scoring PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../common/)
target_link_libraries(scoring common ${EXTERNAL_LIBRARIES})

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/scores.dat DESTINATION .)