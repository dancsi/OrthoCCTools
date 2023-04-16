set(SOURCES 
	BitSet.cpp
)

set(HEADERS
	BB_ColorRSort.h
	BB_GreedyColorSort.h
	BitSet.h
	DegreeAndNumberSort.h
	DegreeSort.h
	GreedyColorSort.h
	KillTimer.h
	MaximumCliqueBase.h
	McrBB.h
	ParallelMaximumClique.h
	SteadyGreedyColorSort.h
	Timer.h
)

find_package(Threads)

add_library(mcqd_para ${SOURCES} ${HEADERS})
target_link_libraries(mcqd_para Threads::Threads)

target_include_directories(mcqd_para PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})