set(COMMON_SOURCES 
	logging.cpp
	PeptideSet.cpp
)

set(COMMON_HEADERS
	logging.h
	SpecialMatrices.h
	MemoryMappedMatrix.h
	PeptideSet.h
	TriangluarIndex.h
	ParallelFor.h
)

add_library(common ${COMMON_SOURCES} ${COMMON_HEADERS})
target_link_libraries(common spdlog::spdlog)

target_include_directories(common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})