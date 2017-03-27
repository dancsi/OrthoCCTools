set(SETBUILDER_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/setbuilder.cpp
)

set(SETBUILDER_ALL ${SETBUILDER_SOURCES} ${SETBUILDER_HEADERS})

add_executable(setbuilder ${SETBUILDER_ALL})

target_link_libraries(setbuilder common scoring flags ${EXTERNAL_LIBRARIES})

install(TARGETS setbuilder RUNTIME DESTINATION .)