add_library(flags INTERFACE )
target_include_directories(flags INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/flags.h)

#target_sources(flags INTERFACE ${SOURCE})