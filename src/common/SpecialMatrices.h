#pragma once

#include <MemoryMappedMatrix.h>

#include <string_view>

template<typename T>
struct SquareMatrix
{
	explicit SquareMatrix(size_t n) : n{n}, ptr { new T[n * n] } {}

	~SquareMatrix() {
		delete ptr;
	}

	T* operator[](size_t index) {
		return ptr + n * index;
	}

	template<typename T>
	static SquareMatrix from_binary(std::string_view path) {
		MemoryMappedMatrix<T> src(path);

		size_t n = src.get_dimensions().first;
		SquareMatrix<T> dst(n);
		
		memcpy(dst[0], src[0], n * n * sizeof(T));

		return dst;
	}

	size_t size() const {
		return n;
	}
private:
	size_t n;
	T* ptr;
};

using score_t = float;
using InteractionMatrix = SquareMatrix<score_t>;
