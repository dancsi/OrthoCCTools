#include "ScoringEnginePotapov.h"

#include "common/experimental_cxx_features.h"

#include <exception>
#include <fstream>
#include <limits>
#include <string>

using std::string_view;

ScoringEnginePotapov::ScoringEnginePotapov(string_view weights_path, int max_peptide_length) : max_peptide_length(max_peptide_length)
{
	using namespace std::experimental::filesystem;
	using namespace std::string_literals;

	path path(weights_path.data());
	if (!exists(path)) {
		auto newpath = current_path() / path.filename();
		if (exists(newpath))
			path = newpath;
		else {
			throw std::runtime_error("Requested path"s + path.string() + "does not exist"s);
		}
	}

	pair_weights.reserve(2486);
	triple_weights.reserve(14093);

	std::ifstream fin(path.string());
	std::string registers, residues;
	float value;
	while (fin >> registers >> residues >> value) {
		if (registers.length() == 2) {
			insert_weight<2>(registers, residues, value, pair_weights, pair_register_map);
		}
		else {
			insert_weight<3>(registers, residues, value, triple_weights, triple_register_map);
		}
	}
	fin.close();

	init_pairs();
	init_triples();
}

template<int k, typename weights_type>
float ScoringEnginePotapov::generic_score(string_view chain1, string_view chain2, std::vector<ResidueTuple<k>>& tuples, weights_type& weight_vec) {
	string_view chains[] = { chain1, chain2 };

	float res = 0;

	auto max_length = std::max(chain1.length(), chain2.length());
	for (const auto& tuple : tuples) {
		if (tuple.max_pos() >= max_length) break;
		auto& weight_array = weight_vec[tuple.weight_array_index];

		char buf[k];
		for (int i = 0; i < k; i++) {
			auto chain_id = tuple.res[i].chain_id;
			auto pos = tuple.res[i].pos;

			if (pos >= chains[chain_id].size() || chains[chain_id][pos] == '-') goto skip;

			buf[i] = chains[chain_id][pos];
		}

		res += weight_array[detail::residues_hash<k>({ buf, k })];
	skip: continue;
	}

	return res;
}

float ScoringEnginePotapov::score(string_view chain1, string_view chain2)
{
	const float w0 = -4.54197f;
	return w0 + generic_score(chain1, chain2, pairs, pair_weights) + generic_score(chain1, chain2, triples, triple_weights);
}

int ScoringEnginePotapov::register_idx(std::string_view registers, std::map<std::string, int>& rmap)
{
	auto it = rmap.find((std::string)registers);
	if (it != rmap.end()) {
		return it->second;
	}

	int next_idx = static_cast<int>(rmap.size());
	rmap.insert({ (std::string)registers, next_idx });

	return next_idx;
}

//we will use this to sort the generated pairs and triples on the maximum position in any chain
template<typename T>
bool max_pos_comparator(T p1, T p2) {
	auto element_key = [](decltype(p1) p) -> decltype(auto) {
		return p.max_pos();
	};

	return element_key(p1) < element_key(p2);
}

void ScoringEnginePotapov::init_pairs()
{
	for (int i = 0; i < max_peptide_length; i++) {
		for (int j = std::max(0, i - 7 + 1); j < std::min(max_peptide_length, i + 7); j++) {
			//ordered registers
			const char ha[] = "FGABCDE";

			//construct pair
			ResiduePair pair{ { (uint16_t)i, 0, ha[i % 7] },{ (uint16_t)j, 1, ha[j % 7] } };

			//order it
			if (i > j) std::swap(pair.res[0], pair.res[1]);

			//if we have some weights for the register sequence, remember them
			std::string pair_registers{ pair.res[0].reg, pair.res[1].reg };
			auto possible_idx = pair_register_map.find(pair_registers);
			if (possible_idx != pair_register_map.end()) {
				pair.weight_array_index = possible_idx->second;
				pairs.push_back(pair);
			}

		}
	}

	std::sort(pairs.begin(), pairs.end(), max_pos_comparator<ResiduePair>);
}

void ScoringEnginePotapov::init_triples()
{
	for (int _i = 0; _i < max_peptide_length; _i++) {
		for (int _j = _i + 1; _j < 2 * max_peptide_length; _j++) {
			for (int _k = std::max(_j + 1, max_peptide_length); _k < 2 * max_peptide_length; _k++) {

				//residue positions in their respective chains
				uint16_t i = _i % max_peptide_length,
					j = _j % max_peptide_length,
					k = _k % max_peptide_length;

				//to which chains do they belong
				uint8_t c_i = _i / max_peptide_length,
					c_j = _j / max_peptide_length,
					c_k = _k / max_peptide_length;

				//if they are more than one heptad apart, there is no interaction
				if (abs(i - j) >= 7 || abs(i - k) >= 7 || abs(k - j) >= 7) continue;

				//if they all belong to the same chain, we don't care
				if (c_i == c_j && c_j == c_k) continue;

				//ordered registers
				const char ha[] = "FGABCDE";

				//form a ResidueTriple
				ResidueTriple triple{ { i, c_i, ha[i % 7] },{ j, c_j, ha[j % 7] },{ k, c_k, ha[k % 7] } };

				//canonicalize triple registers
				std::sort(triple.res.begin(), triple.res.end(), [](auto t1, auto t2) {
					return t1.pos == t2.pos ? t1.chain_id < t2.chain_id : t1.pos < t2.pos;
				});

				//reorder so that pointers from the same chain are adjacent
				if (triple.res[0].pos == triple.res[1].pos &&
					triple.res[0].chain_id == triple.res[2].chain_id) {
					std::swap(triple.res[0], triple.res[1]);
				}

				//we will use this to look up the appropriate weight array
				std::string triad_registers{ triple.res[0].reg, triple.res[1].reg, triple.res[2].reg };

				//only some combinations of chain_ids are allowed
				bool difference1 = triple.res[0].chain_id != triple.res[1].chain_id;
				bool difference2 = triple.res[0].chain_id != triple.res[2].chain_id;

				static const std::map<std::string, std::pair<bool, bool>> allowed_chain_sequences = {
					{ "GDE" ,{ 0,1 } },
					{ "GAE" ,{ 1,1 } },
					{ "DEG" ,{ 1,0 } },
					{ "EGA" ,{ 1,0 } },

					{ "GAD" ,{ 1,0 } },
					{ "ADE" ,{ 1,0 } },
					{ "DGA" ,{ 0,1 } },
					{ "DEA" ,{ 1,1 } },

					{ "AAD" ,{ 1,1 } },
					{ "DDA" ,{ 1,1 } },
				};

				//if we could not find the register sequence, or it is wrong, we don't care
				auto seq = allowed_chain_sequences.find(triad_registers);
				if (seq == allowed_chain_sequences.end() || seq->second != std::make_pair(difference1, difference2)) continue;

				//if we have the corresponding weights, cache that index, and remember the triple
				auto possible_idx = triple_register_map.find(triad_registers);
				if (possible_idx != triple_register_map.end()) {
					triple.weight_array_index = possible_idx->second;
					triples.push_back(triple);
				}
			}
		}
	}

	std::sort(triples.begin(), triples.end(), max_pos_comparator<ResidueTriple>);
}