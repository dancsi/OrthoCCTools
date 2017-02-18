#pragma once

#include <utility>

#include "MemoryMappedMatrix.h"

template<typename T>
struct SquareMatrix : public MemoryMappedMatrix<T> {
public: 
	SquareMatrix(std::string_view path) : MemoryMappedMatrix<T>::MemoryMappedMatrix(path) {}
	SquareMatrix(std::string_view path, size_t n) : MemoryMappedMatrix<T>::MemoryMappedMatrix(path, n, n) {}
};

using score_t = float;
using InteractionMatrix = SquareMatrix<score_t>;