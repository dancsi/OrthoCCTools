set(SOLVER_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/solver.cpp
)

set(SOLVER_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/common.h
	${CMAKE_CURRENT_SOURCE_DIR}/ioutil.h
	${CMAKE_CURRENT_SOURCE_DIR}/options.h
	${CMAKE_CURRENT_SOURCE_DIR}/getopt.h
)

set(SOLVER_ALL ${SOLVER_SOURCES} ${SOLVER_HEADERS})

add_subdirectory(mcqd_para)

add_executable(solver ${SOLVER_ALL})
target_link_libraries(solver common flags MemoryMapped mcqd_para ${EXTERNAL_LIBRARIES})

install(TARGETS solver RUNTIME DESTINATION .)