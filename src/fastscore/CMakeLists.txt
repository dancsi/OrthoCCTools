set(FASTSCORE_SOURCES 
	fastscore.cpp
)

set(FASTSCORE_HEADERS
	options.h
)

add_executable(fastscore ${FASTSCORE_SOURCES} ${FASTSCORE_HEADERS})

target_link_libraries(fastscore common scoring flags MemoryMapped)

if(EMSCRIPTEN)
	target_compile_options(fastscore PRIVATE -pthread)
	target_link_options(fastscore PRIVATE --pre-js ../src/fastscore/wrapper.js -pthread -sPTHREAD_POOL_SIZE=navigator.hardwareConcurrency)
endif()

install(TARGETS fastscore RUNTIME DESTINATION .)