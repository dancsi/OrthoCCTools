set(COMMON_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/logging.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeptideSet.cpp
)

set(COMMON_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/logging.h
	${CMAKE_CURRENT_SOURCE_DIR}/SpecialMatrices.h
	${CMAKE_CURRENT_SOURCE_DIR}/MemoryMappedMatrix.h
	${CMAKE_CURRENT_SOURCE_DIR}/PeptideSet.h
)

add_library(common ${COMMON_SOURCES} ${COMMON_HEADERS})
target_link_libraries(common spdlog::spdlog)

target_include_directories(common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})