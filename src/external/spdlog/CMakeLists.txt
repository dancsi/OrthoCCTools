set(SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/async_logger.h
	${CMAKE_CURRENT_SOURCE_DIR}/common.h
	${CMAKE_CURRENT_SOURCE_DIR}/formatter.h
	${CMAKE_CURRENT_SOURCE_DIR}/logger.h
	${CMAKE_CURRENT_SOURCE_DIR}/spdlog.h
	${CMAKE_CURRENT_SOURCE_DIR}/tweakme.h)

add_library(spdlog INTERFACE )
target_include_directories(spdlog INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../)

#target_sources(spdlog INTERFACE ${SOURCES})