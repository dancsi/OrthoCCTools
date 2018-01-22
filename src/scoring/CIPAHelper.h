#pragma once

#include "ScoringHelper.h"

template<
	typename CIPAImpl
>
struct CIPAHelper {
	using string_view = std::string_view;

	using weights_t = std::array<float, 20 * 20>;
	weights_t c_scores, es_scores;

	void insert_weights(const std::pair<std::string_view, float> weights_to_insert[], size_t length, weights_t& weights) {
		weights.fill(0);
		for (int i = 0; i < length; i++) {
			auto& p = weights_to_insert[i];
			weights[detail::residues_hash<2>(p.first)] = p.second;
		}
	}

	template<typename T, size_t N>
	constexpr size_t array_length(T(&)[N])
	{
		return N;
	}

	void init_c_weights() {
		const auto& weights_to_insert = CIPAImpl::c_weights;
		insert_weights(weights_to_insert, array_length(weights_to_insert), c_scores);
	}
	void init_es_weights() {
		const auto& weights_to_insert = CIPAImpl::es_weights;
		insert_weights(weights_to_insert, array_length(weights_to_insert), es_scores);
	}

	CIPAHelper() {
		init_c_weights();
		init_es_weights();
	}

	float generic_pair_score(const char c1, const char c2, weights_t& weights) {
		char buf[2];
		buf[0] = c1;
		buf[1] = c2;

		return weights[detail::residues_hash<2>({ buf, 2 })];
	}

	float hp_score(const char c) {
		if (!isupper(c)) return 0;
		const float scores[] = { 1.41f, 0.66f, 0.99f, 1.59f, 1.16f, 0.43f, 1.05f, 1.09f, 1.23f, 1.34f, 1.30f, 0.76f, 0.34f, 1.27f, 1.21f, 0.57f, 0.76f, 0.98f, 1.02f, 0.74f };
		return scores[detail::residue_code(c)];
	}
	float c_score(const char c1, const char c2) {
		return generic_pair_score(c1, c2, c_scores);
	}
	float es_score(const char c1, const char c2) {
		return generic_pair_score(c1, c2, es_scores);
	}

	float score(string_view chain1, string_view chain2) {
		auto n1 = chain1.length();
		auto n2 = chain2.length();

		auto n = std::min(n1, n2);

		int n_pairs = 0;
		float hp_sum = 0.f, es_sum = 0.f, c_sum = 0.f;

		for (int i = 0; i < n; i++)
			//fgabcde register
			//0123456
		{
			int reg = i % 7;  //register index (0 = f)

			if (!isupper(chain1[i]) || !isupper(chain2[i])) continue;

			/* helical propensity contribution */
			if ((hp_score(chain1[i])) != 0 && (hp_score(chain2[i]) != 0))  //only count paired residues ('-':'X' pairs contribute nothing)
			{
				hp_sum += hp_score(chain1[i]) + hp_score(chain2[i]);
				n_pairs++;
			}

			/* core contribution */
			if (CIPAImpl::core_position_filter(reg))
			{
				c_sum += c_score(chain1[i], chain2[i]);
			}

			/* electrostatic contribution */
			else if (reg == 1)  //g sites - interact with next e on opposite chain
			{
				if (i + 5 < n2)
				{
					if (!isupper(chain2[i + 5])) continue;
					es_sum += es_score(chain1[i], chain2[i + 5]);
				}
				if (i + 5 < n1)
				{
					if (!isupper(chain1[i + 5])) continue;
					es_sum += es_score(chain2[i], chain1[i + 5]);
				}
			}
		}

		return -CIPAImpl::calculate_score(hp_sum/n_pairs, c_sum, es_sum);
	}
};

