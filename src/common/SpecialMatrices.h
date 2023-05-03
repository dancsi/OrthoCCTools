#pragma once

#include <MemoryMappedMatrix.h>

#include <string_view>

template<typename T>
struct SquareMatrix
{
	explicit SquareMatrix(size_t n) : n{n}, ptr { new T[n * n] } {}

	~SquareMatrix() {
		delete[] ptr;
	}

	T* operator[](size_t index) const {
		return ptr + n * index;
	}

	template<typename U>
	static SquareMatrix from_binary(std::string_view path) {
		MemoryMappedMatrix<U> src(path);

		size_t n = src.get_dimensions().first;
		SquareMatrix<U> dst(n);
		
		memcpy(dst[0], src[0], n * n * sizeof(T));

		return dst;
	}

	void to_binary(std::string_view path) const {
		MemoryMappedMatrix<T> dst(path, n, n);
		memcpy(dst[0], ptr, n * n * sizeof(T));
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
