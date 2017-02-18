#include "ScoringEngine.h"

#include <exception>
#include <fstream>
#include <limits>
#include <string>

using std::string_view;

ScoringEngine::ScoringEngine() {}

float ScoringEngine::score(string_view chain1, string_view chain2)
{
	return 0;
}

ScoringEngine::aligned_score_t ScoringEngine::score(string_view chain1, string_view chain2, alignment_t max_heptad_displacement)
{
	auto left_trimmed = chain2;
	auto right_trimmed = chain2;
	auto len = chain2.length();

	aligned_score_t best_score{ std::numeric_limits<float>::infinity(), 0 };

	for (int displacement = 0; displacement <= max_heptad_displacement; displacement++) {
		aligned_score_t left_trimmed_score = { score(chain1.substr(0, len), left_trimmed), +displacement };
		aligned_score_t right_trimmed_score = { score(chain1.substr(0, len), right_trimmed), -displacement };

		if (left_trimmed_score < best_score) best_score = left_trimmed_score;
		if (right_trimmed_score < best_score) best_score = right_trimmed_score;

		left_trimmed.remove_prefix(7);
		right_trimmed.remove_suffix(7);
		len -= 7;
	}

	return best_score;
}

ScoringEngine::aligned_oriented_score_t ScoringEngine::score(string_view chain1, string_view chain2, alignment_t max_heptad_displacement, Orientation orientation)
{
	static std::string buf;
	aligned_oriented_score_t parallel_score, antiparallel_score;

	if (orientation == Orientation::antiparallel || orientation == Orientation::both) {
		buf.assign(chain2.rbegin(), chain2.rend());
		antiparallel_score = { score(chain1, buf, max_heptad_displacement), Orientation::antiparallel };
	}

	if (orientation == Orientation::parallel || orientation == Orientation::both) {
		parallel_score = { score(chain1, chain2, max_heptad_displacement), Orientation::parallel };
	}

	if (antiparallel_score.score < parallel_score.score) {
		return antiparallel_score;
	}
	else {
		return parallel_score;
	}
}

namespace detail {
	uint8_t residue_code(char r) {
		//The amino acid alphabet is ACDEFGHIKLMNPQRSTVWY
		//We need to map the letters to consecutive numbers
		const uint8_t map[] = { 0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 11, 11, 12, 13, 14, 15, 16, 16, 17, 18, 18, 19, 19 };
		return map[r - 'A'];
	}

	template<>
	int residues_hash<2>(std::string_view residues)
	{
		return 20 * residue_code(residues[0]) + residue_code(residues[1]);
	}

	template<>
	int residues_hash<3>(std::string_view residues)
	{
		return 400 * residue_code(residues[0]) + 20 * residue_code(residues[1]) + residue_code(residues[2]);
	}
}

