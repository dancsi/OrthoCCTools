set(SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/BitSet.cpp
)

set(HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/BB_ColorRSort.h
	${CMAKE_CURRENT_SOURCE_DIR}/BB_GreedyColorSort.h
	${CMAKE_CURRENT_SOURCE_DIR}/BitSet.h
	${CMAKE_CURRENT_SOURCE_DIR}/DegreeAndNumberSort.h
	${CMAKE_CURRENT_SOURCE_DIR}/DegreeSort.h
	${CMAKE_CURRENT_SOURCE_DIR}/GreedyColorSort.h
	${CMAKE_CURRENT_SOURCE_DIR}/KillTimer.h
	${CMAKE_CURRENT_SOURCE_DIR}/MaximumCliqueBase.h
	${CMAKE_CURRENT_SOURCE_DIR}/McrBB.h
	${CMAKE_CURRENT_SOURCE_DIR}/ParallelMaximumClique.h
	${CMAKE_CURRENT_SOURCE_DIR}/SteadyGreedyColorSort.h
	${CMAKE_CURRENT_SOURCE_DIR}/Timer.h
)

find_package(Threads)

add_library(mcqd_para ${SOURCES} ${HEADERS})
target_link_libraries(mcqd_para Threads::Threads)

target_include_directories(mcqd_para PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})