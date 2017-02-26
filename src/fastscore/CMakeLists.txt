set(FASTSCORE_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/fastscore.cpp
)

set(FASTSCORE_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/options.h
)

set(FASTSCORE_ALL ${FASTSCORE_SOURCES} ${FASTSCORE_HEADERS})

add_executable(fastscore ${FASTSCORE_ALL})

target_link_libraries(fastscore common scoring flags MemoryMapped ${EXTERNAL_LIBRARIES})
target_include_directories(fastscore PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../common/)

install(TARGETS fastscore RUNTIME DESTINATION .)