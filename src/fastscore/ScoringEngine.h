#pragma once

#include "common/experimental_cxx_features.h"

#include <array>
#include <map>
#include <tuple>
#include <vector>

struct ScoringEngine {
	const int max_peptide_length;
	ScoringEngine(std::string_view weights_path = "scores.dat", int max_peptide_length = 100);

	enum class Orientation : uint8_t { invalid, parallel, antiparallel, both };
	typedef int8_t alignment_t;

	float score(std::string_view chain1, std::string_view chain2);

	struct aligned_score_t {
		float score;
		alignment_t alignment;

		const bool operator<(aligned_score_t rhs) {
			return score < rhs.score;
		}
	};

	aligned_score_t score(std::string_view chain1, std::string_view chain2, alignment_t max_heptad_displacement);

	struct aligned_oriented_score_t {
		float score;
		alignment_t alignment;
		Orientation orientation;

		aligned_oriented_score_t() : score(std::numeric_limits<float>::infinity()), alignment(0), orientation(Orientation::invalid) {}
		aligned_oriented_score_t(aligned_score_t aligned_score, Orientation orientation) :
			score(aligned_score.score),
			alignment(aligned_score.alignment),
			orientation(orientation)
		{}
	};

	aligned_oriented_score_t score(std::string_view chain1, std::string_view chain2, alignment_t max_heptad_displacement, Orientation orientation);

	std::vector<std::array<float, 20 * 20>> pair_weights;
	std::map<std::string, int> pair_register_map;

	std::vector<std::array<float, 20 * 20 * 20>> triple_weights;
	std::map<std::string, int> triple_register_map;

	int register_idx(std::string_view registers, std::map<std::string, int>& rmap);

	template<int len>
	int residues_hash(std::string_view residues);

	template<size_t len>
	void insert_weight(
		std::string_view registers,
		std::string_view residues,
		float weight,
		std::vector<std::array<float, len>>& weights,
		std::map<std::string, int>& rmap)
	{
		auto idx = register_idx(registers, rmap);
		if (idx >= weights.size()) {
			weights.resize(idx + 1);
		}
		auto h = residues_hash<len>(residues);
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
	float generic_score(std::string_view chain1, std::string_view chain2, std::vector<ResidueTuple<k>>& tuples, weights_type& weight_vec);
};