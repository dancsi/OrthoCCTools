set(COMMON_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/logging.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeptideSet.cpp
)

set(COMMON_HEADERS 
	${CMAKE_CURRENT_SOURCE_DIR}/experimental_cxx_features.h
	${CMAKE_CURRENT_SOURCE_DIR}/logging.h
	${CMAKE_CURRENT_SOURCE_DIR}/SpecialMatrices.h
	${CMAKE_CURRENT_SOURCE_DIR}/MemoryMappedMatrix.h
	${CMAKE_CURRENT_SOURCE_DIR}/PeptideSet.h
)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(common ${COMMON_SOURCES} ${COMMON_HEADERS})
target_link_libraries(common spdlog)
target_compile_features(common PUBLIC cxx_std_17)

target_include_directories(common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})