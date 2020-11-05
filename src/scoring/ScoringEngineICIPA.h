#pragma once

#include "CIPAHelper.h"

#include <array>
#include <string_view>
#include <utility>

using namespace std::literals;

struct ICIPA_core_vert_impl {
	constexpr static std::array<std::pair<std::string_view, float>, 4> c_weights{ {
		{ "NN"sv, -0.290623 }, { "IN"sv, -4.569728 },
		{ "II"sv,  6.792355 }, { "NI"sv, -4.569728 }
	} };

	constexpr static std::array<std::pair<std::string_view, float>, 4> es_weights{ {
		{ "EE"sv, -1.390909 }, { "EK"sv, 4.579879 },
		{ "KK"sv,  0.675038 }, { "KE"sv, 4.579879 }
	} };

	constexpr static std::array<std::pair<std::string_view, float>, 4> cv_weights{ {
		{ "NN"sv, -0.656241 }, { "IN"sv, -2.134678 },
		{ "II"sv,  3.551604 }, { "NI"sv, 3.103324 }
	} };

	constexpr static std::array<std::pair<std::string_view, float>, 0> nterm_c_weights{};

	constexpr static bool core_position_filter(int reg) {
		return (reg == 2);
	}

	constexpr static float calculate_score(CIPAScores scores) {
		return scores.c_sum
			+ scores.es_sum
			+ 1.932004 * scores.num_LL_on_d
			+ -0.274866 * scores.charge_prod
			+ scores.cv_sum;
	}
};

struct ICIPA_nter_core_impl {
	constexpr static std::array<std::pair<std::string_view, float>, 4> c_weights{ {
		{ "NN"sv, -4.077848 }, { "IN"sv, -5.162530 },
		{ "II"sv, 11.497833 }, { "NI"sv, -5.162530 }
	} };

	constexpr static std::array<std::pair<std::string_view, float>, 4> es_weights{ {
		{ "EE"sv, -1.208977 }, { "EK"sv, 4.778290 },
		{ "KK"sv,  0.945598 }, { "KE"sv, 4.778290 }
	} };

	constexpr static std::array<std::pair<std::string_view, float>, 0> cv_weights {};

	constexpr static std::array<std::pair<std::string_view, float>, 4> nterm_c_weights{ {
		{ "NN"sv,  3.843515 }, { "IN"sv, 1.596974 },
		{ "II"sv, -5.440489 }, { "NI"sv, 1.596974 }
	} };

	constexpr static bool core_position_filter(int reg) {
		return (reg == 2);
	}

	constexpr static float calculate_score(CIPAScores scores) {
		return scores.c_sum
			+ scores.es_sum
			+ 2.257456 * scores.num_LL_on_d
			+ -0.282358 * scores.charge_prod
			+ scores.nterm_c;
	}
};

using ScoringEngineICIPACoreVert = CIPAHelper<ICIPA_core_vert_impl>;
using ScoringEngineICIPANterCore = CIPAHelper<ICIPA_nter_core_impl>;
