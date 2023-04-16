set(SOLVER_SOURCES 
	solver.cpp
)

set(SOLVER_HEADERS
	common.h
	ioutil.h
	options.h
	getopt.h
)

add_executable(solver ${SOLVER_SOURCES} ${SOLVER_HEADERS})

add_subdirectory(mcqd_para)
target_link_libraries(solver common flags MemoryMapped mcqd_para)

install(TARGETS solver RUNTIME DESTINATION .)