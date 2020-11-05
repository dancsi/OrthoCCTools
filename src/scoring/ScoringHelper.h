#pragma once

#include <array>
#include <cstring>
#include <map>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace detail {
	inline uint8_t residue_code(const char r) {
		//The amino acid alphabet is ACDEFGHIKLMNPQRSTVWY
		//We need to map the letters to consecutive numbers
		const uint8_t map[] = { 0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 11, 11, 12, 13, 14, 15, 16, 16, 17, 18, 18, 19, 19 };
		return map[r - 'A'];
	}

	template<int len>
	int residues_hash(std::string_view residues);

	template<>
	inline int residues_hash<2>(std::string_view residues)
	{
		return 20 * residue_code(residues[0]) + residue_code(residues[1]);
	}

	template<>
	inline int residues_hash<3>(std::string_view residues)
	{
		return 400 * residue_code(residues[0]) + 20 * residue_code(residues[1]) + residue_code(residues[2]);
	}

	template<typename int_type, int_type x, int k>
	struct pow : std::integral_constant<int_type, x * pow<int_type, x, k - 1>::value> {};

	template<typename int_type, int_type x>
	struct pow<int_type, x, 0> : std::integral_constant<int_type, 1> {};

	template<typename int_type, int_type x, int k>
	constexpr int_type pow_v = pow<int_type, x, k>::value;
}

namespace ScoringOptions {
	enum class Orientation : uint8_t { parallel, antiparallel, both, invalid };
	typedef int32_t alignment_t;

	struct aligned_score_t {
		float score;
		alignment_t alignment;

		const bool operator<(aligned_score_t rhs) {
			return score < rhs.score;
		}
	};

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

	enum class ScoreFunc { potapov, bcipa, qcipa, icipa_core_vert, icipa_nter_core };
}

template<typename ScoringEngine>
struct ScoringHelper {
	using Orientation = ScoringOptions::Orientation;

	using alignment_t = ScoringOptions::alignment_t;
	using aligned_score_t = ScoringOptions::aligned_score_t;
	using aligned_oriented_score_t = ScoringOptions::aligned_oriented_score_t;

	ScoringEngine sc;

	template<typename... Args>
	ScoringHelper(Args&&... args) : sc(std::forward<Args>(args)...) {}

	aligned_score_t score(std::string_view chain1, std::string_view chain2, const std::vector<alignment_t>& alignment, bool truncate)
	{
		const auto max_displacement = 7;
		auto n1 = chain1.length(), n2 = chain2.length();
		auto buffer1_size = n1 + 2 * max_displacement;
		auto buffer2_size = n2 + 2 * max_displacement;

		static thread_local std::vector<char> buf1(buffer1_size, '-'), buf2(buffer2_size, '-');
		if (buf1.size() < buffer1_size) buf1.resize(buffer1_size, '-');
		if (buf2.size() < buffer2_size) buf2.resize(buffer2_size, '-');

		if (truncate) {
			memcpy(buf1.data(), chain1.data(), n1); buf1.resize(n1);
			memcpy(buf2.data(), chain2.data(), n2); buf2.resize(n2);
		}
		else {
			//restore padding
			memset(buf1.data(), '-', max_displacement);
			memset(buf1.data() + max_displacement + n1, '-', max_displacement);
			memset(buf2.data(), '-', max_displacement);
			memset(buf2.data() + max_displacement + n2, '-', max_displacement);

			memcpy(buf1.data() + max_displacement, chain1.data(), n1);
			memcpy(buf2.data() + max_displacement, chain2.data(), n2);
		}
		

		aligned_score_t best_score{ std::numeric_limits<float>::infinity(), 0 };

		const std::string_view original_chain1{ buf1.data(), buf1.size() };
		const std::string_view original_chain2{ buf2.data(), buf2.size() };

		auto shifting_chain1 = original_chain1;
		auto shifting_chain2 = original_chain2;

		auto right_truncated_chain1 = original_chain1;
		auto right_truncated_chain2 = original_chain2;

		alignment_t prev_displacement = 0;

		for (auto displacement : alignment) {
			auto positions_to_trim = displacement - prev_displacement;
			prev_displacement = displacement;

			shifting_chain1.remove_prefix(positions_to_trim);
			shifting_chain2.remove_prefix(positions_to_trim);

			if (truncate) {
				right_truncated_chain1.remove_suffix(positions_to_trim);
				right_truncated_chain2.remove_suffix(positions_to_trim);
			}

			aligned_score_t left_trimmed_score = { sc.score(right_truncated_chain1, shifting_chain2), -displacement };
			aligned_score_t right_trimmed_score = { sc.score(shifting_chain1, right_truncated_chain2), displacement };

			if (left_trimmed_score < best_score) best_score = left_trimmed_score;
			if (right_trimmed_score < best_score) best_score = right_trimmed_score;
		}

		return best_score;
	}

	aligned_oriented_score_t score(std::string_view chain1, std::string_view chain2, const std::vector<alignment_t>& alignment, bool truncate, Orientation orientation)
	{
		static thread_local std::string buf;
		aligned_oriented_score_t parallel_score, antiparallel_score;

		if (orientation == Orientation::antiparallel || orientation == Orientation::both) {
			buf.assign(chain2.rbegin(), chain2.rend());
			antiparallel_score = { score(chain1, buf, alignment, truncate), Orientation::antiparallel };
		}

		if (orientation == Orientation::parallel || orientation == Orientation::both) {
			parallel_score = { score(chain1, chain2, alignment, truncate), Orientation::parallel };
		}

		if (antiparallel_score.score < parallel_score.score) {
			return antiparallel_score;
		}
		else {
			return parallel_score;
		}
	}
};