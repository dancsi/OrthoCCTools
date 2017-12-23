#pragma once

#include "ScoringHelper.h"

struct ScoringEngineQCIPA {
	using string_view = std::string_view;

	using weights_t = std::array<float, 20 * 20>;
	weights_t c_scores, es_scores;

	ScoringEngineQCIPA();

	float score(string_view chain1, string_view chain2);
	float hp_score(const char c);
	float ge_score(char c1, char c2);
};