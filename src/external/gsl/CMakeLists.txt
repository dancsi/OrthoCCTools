add_library(gsl INTERFACE)
target_include_directories(gsl INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/gsl
	${CMAKE_CURRENT_SOURCE_DIR}/gsl_algorithm
	${CMAKE_CURRENT_SOURCE_DIR}/gsl_assert
	${CMAKE_CURRENT_SOURCE_DIR}/gsl_byte
	${CMAKE_CURRENT_SOURCE_DIR}/gsl_util
	${CMAKE_CURRENT_SOURCE_DIR}/multi_span
	${CMAKE_CURRENT_SOURCE_DIR}/span
	${CMAKE_CURRENT_SOURCE_DIR}/string_span)

#target_sources(gsl INTERFACE ${SOURCES})