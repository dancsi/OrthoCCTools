pybind11_add_module(pyccscore pyccscore.cpp)
target_link_libraries(pyccscore PRIVATE scoring)
install(TARGETS pyccscore LIBRARY DESTINATION .)