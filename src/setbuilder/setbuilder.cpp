#include <algorithm>
#include <iostream>
#include <string>

#include "flags.h"

#include "common/PeptideSet.h"

#include "scoring/ScoringEnginePotapov.h"
#include <random>

using namespace std;
namespace fs = std::filesystem;

vector<pair<int, int>> initial_pair_set;
PeptideSet trimmed_set, initial_set;

float initial_scores[1500][1500];

std::string generate_peptide_id()
{
	static int current_id = 0;
	return "P" + std::to_string(++current_id);
}

pair<std::string_view, std::string_view> choose_initial_pair()
{
	static std::mt19937 gen;
	static std::uniform_int_distribution<size_t> dist(0, initial_pair_set.size() / 4);

	int idx = std::min(initial_pair_set.size() - 1, dist(gen));

	auto p = initial_pair_set[idx];
	return { initial_set[p.first].sequence, initial_set[p.second].sequence };
}

pair<std::string_view, std::string_view> choose_next_pair(std::string initial1, std::string initial2, ScoringEnginePotapov& sc)
{
	static std::mt19937 gen;
	static std::bernoulli_distribution coin(0.25);

	initial1.resize(initial1.size() + 7);
	initial2.resize(initial2.size() + 7);

	auto n = trimmed_set.size();
	float best_score = numeric_limits<float>::infinity();
	int best_i = -1, best_j = -1;

	int new_sequence_dest_idx = initial1.size() - 7;
	int scoring_relevant_idx = initial1.size() - 7;

	for (int i = 0; i < n; i++)
	{
		auto& ext1 = trimmed_set[i].sequence;
		std::copy(ext1.begin(), ext1.end(), initial1.begin() + new_sequence_dest_idx);

		for (int j = 0; j < n; j++)
		{
			auto& ext2 = trimmed_set[j].sequence;
			std::copy(ext2.begin(), ext2.end(), initial2.begin() + new_sequence_dest_idx);

			float score = sc.score({ initial1.c_str() + scoring_relevant_idx, 14 }, { initial2.c_str() + scoring_relevant_idx, 14 });
			float relative_score_difference = abs((best_score - score) / score);

			if (score < best_score)
			{
				if (relative_score_difference > 0.5 || coin(gen))
				{
					best_score = score;
					best_i = i;
					best_j = j;
				}
			}
		}
	}
	return { trimmed_set[best_i].sequence, trimmed_set[best_j].sequence };
}

pair<Peptide, Peptide> generate_single(int heptad_count, ScoringEnginePotapov& sc) {
	std::string p1, p2;
	p1.reserve(1 + 7 * heptad_count); p2.reserve(1 + 7 * heptad_count);
	auto p = choose_initial_pair();
	p1.append(p.first); p2.append(p.second);

	for (int h = 1; h < heptad_count; h++)
	{
		auto p = choose_next_pair(p1, p2, sc);
		p1.append(p.first); p2.append(p.second);
	}

	return { {p1, generate_peptide_id(), 'f' - 'a'}, { p2, generate_peptide_id(), 'f' - 'a' } };
}

void print_pair(FILE *fout, std::pair<Peptide, Peptide>& pp, bool flush = false)
{
	fprintf(fout, ">%s\n%s\n", pp.first.id.c_str(), pp.first.sequence.c_str());
	fprintf(fout, ">%s\n%s\n", pp.second.id.c_str(), pp.second.sequence.c_str());
	if (flush)
	{
		fflush(fout);
	}
}

int main(int argc, char **argv) {
	const flags::args args(argc, argv);
	auto positional = args.positional();

	string input_path(positional[0].data());
	int count = stoi(positional[1].data());
	int heptad_count = stoi(positional[2].data());

	initial_set.read(input_path);

	for (auto& p : initial_set) {
		trimmed_set.push_back(p.remove_padding());
	}

	ScoringEnginePotapov sc{};

	cout << "Initial scoring\n";

	for (int i = 0; i < trimmed_set.size(); i++) {
		for (int j = i; j < trimmed_set.size(); j++) {
			initial_pair_set.push_back(make_pair(i, j));
			initial_pair_set.push_back(make_pair(j, i));
			initial_scores[i][j] = initial_scores[j][i] = sc.score(initial_set[i].sequence, initial_set[j].sequence);
		}
	}

	cout << "Initial scoring done\n";

	sort(initial_pair_set.begin(), initial_pair_set.end(), [](auto p1, auto p2) {
		return initial_scores[p1.first][p1.second] < initial_scores[p2.first][p2.second];
	});

	vector<pair<Peptide, Peptide>> results;

	FILE *fout = fopen("chained.fasta", "w");
	for (int i = 0; i < count; i++) {
		auto p = generate_single(heptad_count, sc);
		results.push_back(p);
		print_pair(stdout, p, true);
		print_pair(fout, p, true);
	}
	fclose(fout);

	return 0;
}