#pragma once

#include "common/experimental_cxx_features.h"

#include <array>
#include <map>
#include <tuple>
#include <vector>

namespace detail {
	uint8_t residue_code(char r);

	template<int len>
	int residues_hash(std::string_view residues);

	template<typename int_type, int_type x, int k>
	struct pow : std::integral_constant<int_type, x * pow<int_type, x, k-1>::value> {};

	template<typename int_type, int_type x>
	struct pow<int_type, x, 0> : std::integral_constant<int_type, 1> {};

	template<typename int_type, int_type x, int k>
	constexpr int_type pow_v = pow<int_type, x, k>::value;
}

struct ScoringEngine {
	ScoringEngine();

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

	template<int len>
	int residues_hash(std::string_view residues) {
		return detail::residues_hash<len>(residues);
	}

};