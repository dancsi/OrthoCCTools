#include "ScoringEngineBCIPA.h"

#include <cctype>

ScoringEngineBCIPA::ScoringEngineBCIPA() {
	init_c_weights();
	init_es_weights();
}

template<typename T, size_t N>
size_t array_length(T(&)[N])
{
	return N;
}

void ScoringEngineBCIPA::init_c_weights() {
	std::pair<std::string, float> weights_to_insert[] = { { "II", -1.5 }, { "LL", -1.5 }, { "VI", -1.5 }, { "IV", -1.5 }, { "VV", -1.0 }, { "VL", -1.0 }, { "IL", -1.0 }, { "IR", -1.0 }, { "IK", -1.0 }, { "LV", -1.0 }, { "LI", -1.0 }, { "RI", -1.0 }, { "KI", -1.0 }, { "IA", -0.5 }, { "LA", -0.5 }, { "VA", -0.5 }, { "NN", -0.5 }, { "IN", -0.5 }, { "IT", -0.5 }, { "LK", -0.5 }, { "LT", -0.5 }, { "RR", -0.5 }, { "AI", -0.5 }, { "AL", -0.5 }, { "AV", -0.5 }, { "NI", -0.5 }, { "TI", -0.5 }, { "KL", -0.5 }, { "TL", -0.5 }, { "VT", +0.5 }, { "TV", +0.5 } };
	insert_weights(weights_to_insert, array_length(weights_to_insert), c_scores);
}

void ScoringEngineBCIPA::init_es_weights() {
	std::pair<std::string, float> weights_to_insert[] = { {"RE", -2.0}, { "ER", -2.0 }, { "KE", -1.5 }, { "KQ", -1.5 }, { "RQ", -1.5 }, { "QQ", -1.5 }, { "EK", -1.5 }, { "QK", -1.5 }, { "QR", -1.5 }, { "QE", -1.0 }, { "EQ", -1.0 }, { "QA", -0.5 }, { "RA", -0.5 }, { "KD", -0.5 }, { "RD", -0.5 }, { "KL", -0.5 }, { "TL", -0.5 }, { "RK", -0.5 }, { "AQ", -0.5 }, { "AR", -0.5 }, { "DK", -0.5 }, { "DR", -0.5 }, { "LK", -0.5 }, { "LT", -0.5 }, { "KR", -0.5 }, { "EE", +0.5 }, { "KK", +0.5 }, { "RR", +0.5 }, { "DD", +1.0 }, { "DE", +1.0 }, { "TR", +1.0 }, { "ED", +1.0 }, { "RT", +1.0 } };
	insert_weights(weights_to_insert, array_length(weights_to_insert), es_scores);
}

void ScoringEngineBCIPA::insert_weights(std::pair<std::string, float> weights_to_insert[], size_t length, weights_t& weights) {
	for (int i = 0; i < length; i++) {
		auto& p = weights_to_insert[i];
		weights[detail::residues_hash<2>(p.first)] = p.second;
	}
}

float ScoringEngineBCIPA::score(string_view chain1, string_view chain2) /*
 * assumes sequences are aligned, starting with f position
 * assumes uppercase characters (case sensitive!)
 * all non-aa characters (including whitespace) assumed to take up space
   in alignment, but contribute nothing to the score
 */
{
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
		if ((reg == 2) || (reg == 5))  //a,d sites
		{
			c_sum += c_score(chain1[i], chain2[i]);
		}

		/* electrostatic contribution */
		else if (reg == 1)  //g sites - interact with next e on opposite chain
		{
			if (i + 5 < n2)
			{
				es_sum += es_score(chain1[i], chain2[i + 5]);
			}
			if (i + 5 < n1)
			{
				es_sum += es_score(chain2[i], chain1[i + 5]);
			}
		}
	}

	return 81.3256f * hp_sum / n_pairs - 10.5716f * c_sum - 4.7771f * es_sum - 29.1320f - 273;
}

float ScoringEngineBCIPA::hp_score(const char c) {
	if (!isupper(c)) return 0;
	const float scores[] = { 1.41f, 0.66f, 0.99f, 1.59f, 1.16f, 0.43f, 1.05f, 1.09f, 1.23f, 1.34f, 1.30f, 0.76f, 0.34f, 1.27f, 1.21f, 0.57f, 0.76f, 0.98f, 1.02f, 0.74f };
	return scores[detail::residue_code(c)];
}

float ScoringEngineBCIPA::c_score(const char c1, const char c2) {
	return generic_pair_score(c1, c2, c_scores);
}

float ScoringEngineBCIPA::es_score(const char c1, const char c2) {
	return generic_pair_score(c1, c2, es_scores);
}

inline float ScoringEngineBCIPA::generic_pair_score(const char c1, const char c2, weights_t& weights) {
	static char buf[2];
	buf[0] = c1;
	buf[1] = c2;

	return weights[detail::residues_hash<2>({ buf, 2 })];
}