set(SETBUILDER_SOURCES 
	setbuilder.cpp
)

add_executable(setbuilder ${SETBUILDER_SOURCES})

target_link_libraries(setbuilder common scoring flags)

install(TARGETS setbuilder RUNTIME DESTINATION .)