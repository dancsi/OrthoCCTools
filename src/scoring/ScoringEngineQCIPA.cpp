#include "ScoringEngineQCIPA.h"

#include <algorithm>
#include <cctype>

ScoringEngineQCIPA::ScoringEngineQCIPA() {
	init_c_weights();
	init_es_weights();
}

template<typename T, size_t N>
inline size_t array_length(T(&)[N])
{
	return N;
}

void ScoringEngineQCIPA::init_c_weights() {
	std::pair<std::string, float> weights_to_insert[] = { {"II", -1.75}, {"IN", 11.78}, {"NI", 11.78}, {"NN", -5.24} };
	insert_weights(weights_to_insert, array_length(weights_to_insert), c_scores);
}

void ScoringEngineQCIPA::init_es_weights() {
	std::pair<std::string, float> weights_to_insert[] = { {"EE", -11.3}, {"EK", -0.97}, {"KE", -0.97}, {"KK", -76.22} };
	insert_weights(weights_to_insert, array_length(weights_to_insert), es_scores);
}

void ScoringEngineQCIPA::insert_weights(std::pair<std::string, float> weights_to_insert[], size_t length, weights_t& weights) {
	weights.fill(0);
	for (int i = 0; i < length; i++) {
		auto& p = weights_to_insert[i];
		weights[detail::residues_hash<2>(p.first)] = p.second;
	}
}

float ScoringEngineQCIPA::score(string_view chain1, string_view chain2) /*
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
		if (reg == 2)  //a,a' sites
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

	return -(4.16 * hp_sum / n_pairs + c_sum + es_sum + 30.18);
}

float ScoringEngineQCIPA::hp_score(const char c) {
	if (!isupper(c)) return 0;
	const float scores[] = { 1.41f, 0.66f, 0.99f, 1.59f, 1.16f, 0.43f, 1.05f, 1.09f, 1.23f, 1.34f, 1.30f, 0.76f, 0.34f, 1.27f, 1.21f, 0.57f, 0.76f, 0.98f, 1.02f, 0.74f };
	return scores[detail::residue_code(c)];
}

float ScoringEngineQCIPA::c_score(const char c1, const char c2) {
	return generic_pair_score(c1, c2, c_scores);
}

float ScoringEngineQCIPA::es_score(const char c1, const char c2) {
	return generic_pair_score(c1, c2, es_scores);
}

inline float ScoringEngineQCIPA::generic_pair_score(const char c1, const char c2, weights_t& weights) {
	char buf[2];
	buf[0] = c1;
	buf[1] = c2;

	return weights[detail::residues_hash<2>({ buf, 2 })];
}
