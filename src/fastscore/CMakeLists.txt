set(FASTSCORE_SOURCES 
	fastscore.cpp
)

set(FASTSCORE_HEADERS
	options.h
)

add_executable(fastscore ${FASTSCORE_SOURCES} ${FASTSCORE_HEADERS})

target_link_libraries(fastscore common scoring flags MemoryMapped)

if(EMSCRIPTEN)
	target_link_options(fastscore PRIVATE -lnoderawfs.js -s NODERAWFS=1)
endif()

install(TARGETS fastscore RUNTIME DESTINATION .)