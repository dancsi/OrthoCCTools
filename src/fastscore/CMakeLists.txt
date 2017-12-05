set(FASTSCORE_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/fastscore.cpp
)

set(FASTSCORE_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/options.h
)

set(FASTSCORE_ALL ${FASTSCORE_SOURCES} ${FASTSCORE_HEADERS})

add_executable(fastscore ${FASTSCORE_ALL})

find_package(OpenMP REQUIRED)
target_compile_options(fastscore PRIVATE ${OpenMP_CXX_FLAGS})

target_link_libraries(fastscore common scoring flags MemoryMapped ${EXTERNAL_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
target_include_directories(fastscore PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../common/)

install(TARGETS fastscore RUNTIME DESTINATION .)

add_custom_command(TARGET fastscore POST_BUILD
				   COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/../scoring/scores.dat" $<TARGET_FILE_DIR:fastscore>)