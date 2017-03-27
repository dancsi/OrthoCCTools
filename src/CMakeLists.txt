INCLUDE (CheckIncludeFileCXX)

macro(CHECK_EXPERIMENTAL_INCLUDE INCLUDE)
	CHECK_INCLUDE_FILE_CXX("${INCLUDE}" "HAVE_${INCLUDE}")
	CHECK_INCLUDE_FILE_CXX("experimental/${INCLUDE}" "HAVE_EXPERIMENTAL_${INCLUDE}")
endmacro()

CHECK_EXPERIMENTAL_INCLUDE("optional")
CHECK_EXPERIMENTAL_INCLUDE("string_view")
CHECK_EXPERIMENTAL_INCLUDE("filesystem")
CHECK_INCLUDE_FILE_CXX("experimental/array" "HAVE_EXPERIMENTAL_array")

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(external)
add_subdirectory(common)
add_subdirectory(scoring)
add_subdirectory(fastscore)
add_subdirectory(solver)
add_subdirectory(setbuilder)