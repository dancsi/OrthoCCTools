#pragma once

#include <algorithm>
#include <map>
#include <vector>

#include "ScoringHelper.h"

struct ScoringEnginePotapov {
	using string_view = std::string_view;

public:
	const int max_peptide_length;

	ScoringEnginePotapov(string_view weights_path = "scores.dat", int max_peptide_length = 100);

	float score(string_view chain1, string_view chain2);

private:

	template<size_t k>
	void insert_weight(
		string_view registers,
		string_view residues,
		float weight,
		std::vector<std::array<float, detail::pow<size_t, 20, k>::value>>& weights,
		std::map<std::string, int>& rmap)
	{
		auto idx = register_idx(registers, rmap);
		if (idx >= weights.size()) {
			weights.resize(idx + 1);
		}
		auto h = detail::residues_hash<k>(residues);
		weights[idx][h] = weight;
	}

	struct ResiduePointer {
		uint16_t pos;
		uint8_t chain_id;
		char reg;
	};

	template<int len>
	struct ResidueTuple {
		std::array<ResiduePointer, len> res;
		int weight_array_index;
		uint16_t max_pos_cache;

		inline uint16_t max_pos() const {
			return max_pos_cache;
		}

		ResidueTuple(std::initializer_list<ResiduePointer> list) : max_pos_cache(0) {
			std::move(list.begin(), list.end(), res.begin());
			max_pos_cache = std::max_element(res.begin(), res.end(),
				[](auto t1, auto t2) {
				return t1.pos < t2.pos;
			})->pos;
		}
	};

	typedef ResidueTuple<2> ResiduePair;
	typedef ResidueTuple<3> ResidueTriple;

	std::vector<ResiduePair> pairs;
	std::vector<ResidueTriple> triples;

	void init_pairs();
	void init_triples();

	template<int k, typename weights_type>
	float generic_score(string_view chain1, string_view chain2, std::vector<ResidueTuple<k>>& tuples, weights_type& weight_vec);

	std::vector<std::array<float, 20 * 20>> pair_weights;
	std::map<std::string, int> pair_register_map;

	std::vector<std::array<float, 20 * 20 * 20>> triple_weights;
	std::map<std::string, int> triple_register_map;

	int register_idx(string_view registers, std::map<std::string, int>& rmap);
};