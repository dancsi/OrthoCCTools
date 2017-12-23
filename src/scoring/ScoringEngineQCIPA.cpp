#include "ScoringEngineQCIPA.h"

#include <algorithm>
#include <cctype>

ScoringEngineQCIPA::ScoringEngineQCIPA() {
}

template<typename T, size_t N>
inline size_t array_length(T(&)[N])
{
	return N;
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
	float hp_sum = 0, a_sum = 0, ge_sum;

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
		if (reg == 2)  //sites a, a'
		{
			char c1 = chain1[i];
			char c2 = chain2[i];
			if (c1 > c2) std::swap(c1, c2);
			if (c1 == 'I') {
				if (c2 == 'I') {
					a_sum += -1.75;
				}
				else if (c2 == 'N') {
					a_sum += 11.78;
				}
			}
			else if (c1 == 'N' && c2 == 'N') {
				a_sum += -5.24;
			}
		}

		/* electrostatic contribution */
		else if (reg == 1)  //g sites - interact with next e on opposite chain
		{
			if (i + 5 < n2)
			{
				if (!isupper(chain2[i + 5])) continue;
				ge_sum += ge_score(chain1[i], chain2[i + 5]);
			}
			if (i + 5 < n1)
			{
				if (!isupper(chain1[i + 5])) continue;
				ge_sum += ge_score(chain2[i], chain1[i + 5]);
			}
		}
	}

	return -(4.16 * hp_sum / n_pairs + a_sum + ge_sum + 30.18);
}

float ScoringEngineQCIPA::hp_score(const char c) {
	if (!isupper(c)) return 0;
	const float scores[] = { 1.41f, 0.66f, 0.99f, 1.59f, 1.16f, 0.43f, 1.05f, 1.09f, 1.23f, 1.34f, 1.30f, 0.76f, 0.34f, 1.27f, 1.21f, 0.57f, 0.76f, 0.98f, 1.02f, 0.74f };
	return scores[detail::residue_code(c)];
}

float ScoringEngineQCIPA::ge_score(char c1, char c2)
{
	float ret = 0;
	if (c1 > c2) std::swap(c1, c2);
	if (c1 == 'E') {
		if (c2 == 'E') {
			ret += -11.3;
		}
		else if (c2 == 'K') {
			ret += -0.97;
		}
	}
	else if (c1 == 'K' && c2 == 'K') {
		ret += -76.22;
	}

	return ret;
}
