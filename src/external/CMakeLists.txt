add_subdirectory(flags)
add_subdirectory(MemoryMapped)
add_subdirectory(spdlog)

if(NOT EMSCRIPTEN)
	add_subdirectory(pybind11)
endif()