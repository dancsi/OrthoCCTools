#pragma once

#include <string>
#include <utility>
#include <ostream>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include "MemoryMapped.h"

struct FileHeader
{
	uint16_t version;
	uint16_t flags;
	uint32_t element_size;
	uint64_t n, m;
	uint64_t offset;
};

template <typename ValueType>
class MemoryMappedMatrix
{
private:
	size_t n, m, offset;
	MemoryMapped file;

	FileHeader* getHeaderPointer()
	{
		return reinterpret_cast<FileHeader*> (file.getData());
	}
	ValueType* getDataPointer()
	{
		return reinterpret_cast<ValueType*> (file.getData() + offset);
	}

public:

	MemoryMappedMatrix(std::string_view path)
	{
		file.open(path, MemoryMapped::WholeFile);
		FileHeader *ptr = getHeaderPointer();
		assert(sizeof(ValueType) == ptr->element_size);
		n = ptr->n;
		m = ptr->m;
		offset = ptr->offset;
	}

	MemoryMappedMatrix(std::string_view path, size_t n, size_t m) :n(n), m(m), offset(sizeof(FileHeader))
	{
		size_t fsize = offset + n * m * sizeof(ValueType);

		file.open(path, fsize);

		FileHeader *ptr = getHeaderPointer();

		{
			ptr->version = 1;
			ptr->flags = 0;
			ptr->element_size = sizeof(ValueType);
			ptr->n = n;
			ptr->m = m;
			ptr->offset = offset;
		}
	}

	~MemoryMappedMatrix()
	{
		file.close();
	}

	std::pair<size_t, size_t> get_dimensions()
	{
		return std::make_pair(n, m);
	}

	friend std::ostream& operator<<(std::ostream& out, MemoryMappedMatrix& mat)
	{
		return out << "MemoryMappedMatrix{ n=" << mat.n << ", m=" << mat.m << " }";
	}

	ValueType *operator[](size_t index)
	{
		return getDataPointer() + index * m;
	}
};


