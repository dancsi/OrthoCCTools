#pragma once

#include <string_view>

struct PotapovScore
{
	std::string_view registers;
	std::string_view residues;
	float value;
};

extern PotapovScore potapov_scores[16579];