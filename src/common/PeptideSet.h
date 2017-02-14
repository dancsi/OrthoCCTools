#pragma once

#include "experimental_cxx_features.h"

#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

class Peptide
{
public:
	const std::string sequence, id;
	Peptide(std::string& sequence, std::string& id) : sequence(sequence), id(id)
	{
		if (id.empty())
			throw std::invalid_argument("peptide id must not be empty");
	}

};

class PeptideSet
{
private:
	std::vector<Peptide> storage;

public:
	PeptideSet() {}
	PeptideSet(std::string_view path)
	{
		read(path);
	}

	void read(std::string_view path);
	void write(std::string& path);

	Peptide& operator[] (size_t index)
	{
		return storage[index];
	}

	size_t size()
	{
		return storage.size();
	}

	decltype(storage)::iterator begin()
	{
		return storage.begin();
	}

	decltype(storage)::iterator end()
	{
		return storage.end();
	}
};