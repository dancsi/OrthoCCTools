#pragma once

#include "ScoringHelper.h"

struct ScoringEngineBCIPA {
	using string_view = std::string_view;

	using weights_t = std::array<float, 20 * 20>;
	weights_t c_scores, es_scores;

	ScoringEngineBCIPA();
	void init_c_weights();
	void init_es_weights();
	void insert_weights(std::pair<std::string, float> weights_to_insert[], size_t length, weights_t& weights);

	float score(string_view chain1, string_view chain2);

	float hp_score(const char c);
	float c_score(const char c1, const char c2);
	float es_score(const char c1, const char c2);

	float generic_pair_score(const char c1, const char c2, weights_t& weights);
};