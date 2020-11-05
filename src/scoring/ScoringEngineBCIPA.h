#pragma once

#include "CIPAHelper.h"

#include <array>
#include <string_view>
#include <utility>

using namespace std::literals;

struct BCIPAImpl {
	constexpr static std::array<std::pair<std::string_view, float>, 31>  c_weights { {
		{ "II"sv, -1.5 },{ "LL"sv, -1.5 },{ "VI"sv, -1.5 },{ "IV"sv, -1.5 },
		{ "VV"sv, -1.0 },{ "VL"sv, -1.0 },{ "IL"sv, -1.0 },{ "IR"sv, -1.0 },
		{ "IK"sv, -1.0 },{ "LV"sv, -1.0 },{ "LI"sv, -1.0 },{ "RI"sv, -1.0 },
		{ "KI"sv, -1.0 },{ "IA"sv, -0.5 },{ "LA"sv, -0.5 },{ "VA"sv, -0.5 },
		{ "NN"sv, -0.5 },{ "IN"sv, -0.5 },{ "IT"sv, -0.5 },{ "LK"sv, -0.5 },
		{ "LT"sv, -0.5 },{ "RR"sv, -0.5 },{ "AI"sv, -0.5 },{ "AL"sv, -0.5 },
		{ "AV"sv, -0.5 },{ "NI"sv, -0.5 },{ "TI"sv, -0.5 },{ "KL"sv, -0.5 },
		{ "TL"sv, -0.5 },{ "VT"sv, +0.5 },{ "TV"sv, +0.5 } 
	} };

	constexpr static std::array<std::pair<std::string_view, float>, 33> es_weights { {
		{ "RE"sv, -2.0 },{ "ER"sv, -2.0 },{ "KE"sv, -1.5 },{ "KQ"sv, -1.5 },
		{ "RQ"sv, -1.5 },{ "QQ"sv, -1.5 },{ "EK"sv, -1.5 },{ "QK"sv, -1.5 },
		{ "QR"sv, -1.5 },{ "QE"sv, -1.0 },{ "EQ"sv, -1.0 },{ "QA"sv, -0.5 },
		{ "RA"sv, -0.5 },{ "KD"sv, -0.5 },{ "RD"sv, -0.5 },{ "KL"sv, -0.5 },
		{ "TL"sv, -0.5 },{ "RK"sv, -0.5 },{ "AQ"sv, -0.5 },{ "AR"sv, -0.5 },
		{ "DK"sv, -0.5 },{ "DR"sv, -0.5 },{ "LK"sv, -0.5 },{ "LT"sv, -0.5 },
		{ "KR"sv, -0.5 },{ "EE"sv, +0.5 },{ "KK"sv, +0.5 },{ "RR"sv, +0.5 },
		{ "DD"sv, +1.0 },{ "DE"sv, +1.0 },{ "TR"sv, +1.0 },{ "ED"sv, +1.0 },
		{ "RT"sv, +1.0 } 
	} };

	constexpr static std::array<std::pair<std::string_view, float>, 0> cv_weights{};

	constexpr static std::array<std::pair<std::string_view, float>, 0> nterm_c_weights{};

	constexpr static bool core_position_filter(int reg) {
		return (reg == 2) || (reg == 5);
	}

	constexpr static float calculate_score(CIPAScores scores) {
		return 
			81.3256f * scores.avg_hp_sum 
			- 10.5716f * scores.c_sum 
			- 4.7771f * scores.es_sum 
			- 29.1320f - 273;
	}
};

using ScoringEngineBCIPA = CIPAHelper<BCIPAImpl>;