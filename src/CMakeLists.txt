INCLUDE (CheckIncludeFileCXX)

macro(CHECK_EXPERIMENTAL_INCLUDE INCLUDE)
	CHECK_INCLUDE_FILE_CXX("${INCLUDE}" "HAVE_${INCLUDE}")
	CHECK_INCLUDE_FILE_CXX("experimental/${INCLUDE}" "HAVE_EXPERIMENTAL_${INCLUDE}")
endmacro()

CHECK_EXPERIMENTAL_INCLUDE("optional")
CHECK_EXPERIMENTAL_INCLUDE("string_view")
CHECK_EXPERIMENTAL_INCLUDE("filesystem")
CHECK_INCLUDE_FILE_CXX("experimental/array" "HAVE_EXPERIMENTAL_array")

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/common/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/common/config.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(external)
add_subdirectory(common)
add_subdirectory(fastscore)