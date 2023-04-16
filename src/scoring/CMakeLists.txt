set(SCORING_SOURCES 
	ScoringEnginePotapov.cpp
	PotapovScores.cpp
)

set(SCORING_HEADERS
	CIPAHelper.h
	ScoringEngineBCIPA.h
	ScoringEngineQCIPA.h
	ScoringEngineICIPA.h
	ScoringEnginePotapov.h
	ScoringHelper.h
	PotapovScores.h
)

add_library(scoring ${SCORING_SOURCES} ${SCORING_HEADERS})

target_link_libraries(scoring common)
set_target_properties(scoring PROPERTIES POSITION_INDEPENDENT_CODE ON)