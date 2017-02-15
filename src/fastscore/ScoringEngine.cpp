#include "ScoringEngine.h"

#include <exception>
#include <fstream>
#include <limits>
#include <string>

using std::string_view;

ScoringEngine::ScoringEngine(string_view weights_path, int max_peptide_length) : max_peptide_length(max_peptide_length)
{
	using namespace std::experimental::filesystem;
	using namespace std::string_literals;

	path path(weights_path.data());
	if (!exists(path)) {
		throw std::runtime_error("Requested path"s + path.string() + "does not exist"s);
	}

	pair_weights.reserve(2486);
	triple_weights.reserve(14093);

	std::ifstream fin(path.string());
	std::string registers, residues;
	float value;
	while (fin >> registers >> residues >> value) {
		if (registers.length() == 2) {
			insert_weight(registers, residues, value, pair_weights, pair_register_map);
		}
		else {
			insert_weight(registers, residues, value, triple_weights, triple_register_map);
		}
	}
	fin.close();

	init_pairs();
	init_triples();
}

template<int k, typename weights_type>
float ScoringEngine::generic_score(string_view chain1, string_view chain2, std::vector<ResidueTuple<k>>& tuples, weights_type& weight_vec) {
	string_view chains[] = { chain1, chain2 };

	float res = 0;

	auto max_length = std::max(chain1.length(), chain2.length());
	for (const auto& tuple : tuples) {
		if (tuple.max_pos() >= max_length) break;
		auto& weight_array = weight_vec[tuple.weight_array_index];

		static char buf[k];
		for (int i = 0; i < k; i++) {
			auto chain_id = tuple.res[i].chain_id;
			auto pos = tuple.res[i].pos;

			if (pos >= chains[chain_id].size() || chains[chain_id][pos] == '-') goto skip;

			buf[i] = chains[chain_id][pos];
		}

		res += weight_array[residues_hash<k>({ buf, k })];
	skip: continue;
	}

	return res;
}

float ScoringEngine::score(string_view chain1, string_view chain2)
{
	const float w0 = -4.54197f;
	return w0 + generic_score(chain1, chain2, pairs, pair_weights) + generic_score(chain1, chain2, triples, triple_weights);
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

int ScoringEngine::register_idx(std::string_view registers, std::map<std::string, int>& rmap)
{
	auto it = rmap.find((std::string)registers);
	if (it != rmap.end()) {
		return it->second;
	}

	int next_idx = rmap.size();
	rmap.insert({ (std::string)registers, next_idx });

	return next_idx;
}

uint8_t residue_code(char r) {
	//The amino acid alphabet is ACDEFGHIKLMNPQRSTVWY
	//We need to map the letters to consecutive numbers
	static const uint8_t map[] = { 0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 11, 11, 12, 13, 14, 15, 16, 16, 17, 18, 18, 19, 19 };
	return map[r - 'A'];
}

template<int len>
int ScoringEngine::residues_hash(std::string_view residues)
{
	int r = 0;
	for (auto c : residues) {
		r = 20 * r + residue_code(c);
	}
	return r;
}

template<>
int ScoringEngine::residues_hash<2>(std::string_view residues)
{
	return 20 * residue_code(residues[0]) + residue_code(residues[1]);
}

template<>
int ScoringEngine::residues_hash<3>(std::string_view residues)
{
	return 400 * residue_code(residues[0]) + 20 * residue_code(residues[1]) + residue_code(residues[2]);
}

//we will use this to sort the generated pairs and triples on the maximum position in any chain
template<typename T>
bool max_pos_comparator(T p1, T p2) {
	auto element_key = [](decltype(p1) p) -> decltype(auto) {
		return p.max_pos();
	};

	return element_key(p1) < element_key(p2);
}

void ScoringEngine::init_pairs()
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

void ScoringEngine::init_triples()
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
				ResidueTriple triple{ { i, c_i, ha[i % 7] }, { j, c_j, ha[j % 7] }, { k, c_k, ha[k % 7] } };

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
					{"GDE" ,{0,1}},
					{"GAE" ,{1,1}},
					{"DEG" ,{1,0}},
					{"EGA" ,{1,0}},

					{"GAD" ,{1,0}},
					{"ADE" ,{1,0}},
					{"DGA" ,{0,1}},
					{"DEA" ,{1,1}},

					{"AAD" ,{1,1}},
					{"DDA" ,{1,1}},
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
