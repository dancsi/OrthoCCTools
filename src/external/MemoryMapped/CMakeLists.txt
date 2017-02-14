set(SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/MemoryMapped.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MemoryMapped.h
)

if(WIN32)
	set(SOURCES ${SOURCES} mman.h mman.c)
endif()

add_library(MemoryMapped ${SOURCES})
target_include_directories(MemoryMapped SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})