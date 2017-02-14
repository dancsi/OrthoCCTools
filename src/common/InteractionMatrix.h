#pragma once

#include <utility>

#include "MemoryMappedMatrix.h"

using score_t = float;

class InteractionMatrix: public MemoryMappedMatrix<score_t>
{
public:
	InteractionMatrix(std::string_view path) : MemoryMappedMatrix(path) {}
	InteractionMatrix(std::string_view path, size_t n_peptides) : MemoryMappedMatrix(path, n_peptides, n_peptides) {}
};
