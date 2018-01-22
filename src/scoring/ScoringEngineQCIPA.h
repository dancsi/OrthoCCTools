#pragma once

#include "CIPAHelper.h"

#include <array>
#include <string_view>
#include <utility>

using namespace std::literals;

struct QCIPAImpl {
	constexpr static std::pair<std::string_view, float> c_weights[] = {
		{ "II"sv, -1.75 },{ "IN"sv, 11.78 },{ "NI"sv, 11.78 },{ "NN"sv, -5.24 }
	};

	constexpr static std::pair<std::string_view, float> es_weights[] = {
		{ "EE"sv, -11.3 },{ "EK"sv, -0.97 },{ "KE"sv, -0.97 },{ "KK"sv, -76.22 } 
	};

	constexpr static bool core_position_filter(int reg) {
		return (reg == 2);
	}

	constexpr static float calculate_score(float avg_hp_sum, float c_sum, float es_sum) {
		return (4.16 * avg_hp_sum + c_sum + es_sum + 30.18);
	}
};

using ScoringEngineQCIPA = CIPAHelper<QCIPAImpl>;