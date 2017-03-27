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
	const short starting_register;

	Peptide(std::string_view sequence, std::string_view id, short starting_register) : sequence(sequence), id(id), starting_register(starting_register)
	{
		if (id.empty())
			throw std::invalid_argument("peptide id must not be empty");
	}

	Peptide remove_padding() {
		size_t n = sequence.length();
		size_t leading_padding = 0, trailing_padding = 0;

		while (leading_padding < n && sequence[leading_padding] == '-') leading_padding++;
		while (trailing_padding < n && sequence[n - trailing_padding - 1] == '-') trailing_padding++;

		return Peptide(
			sequence.substr(leading_padding, n - (leading_padding + trailing_padding)),
			id + '\'',
			(starting_register + leading_padding) % 7
		);
	}
};

class PeptideSet : public std::vector<Peptide>
{

public:
	PeptideSet() {}
	PeptideSet(std::experimental::filesystem::path path)
	{
		read(path);
	}

	void read(std::experimental::filesystem::path path);
	void write(std::string_view path);
};