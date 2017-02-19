#pragma once

#include "flags.h"
#include "ScoringHelper.h"

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace fs = std::experimental::filesystem;

struct Options {
	const flags::args args;

	void print_usage_and_exit() {
		puts(R"(
USAGE: fastscore INPUT [OPTIONS...]
Available options:
    --basename=PATH                 specify output base name
    --max-heptad-displacement=NUM   try shifting the peptides left and right by up to this many heptads
    --alignment=LIST                try these alignments
    --orientation={parallel, antiparallel, both}
    --score-func={potapov, bcipa}   choose scoring function
)");
		exit(1);
	}

	fs::path fasta_path;
	std::string basename;
	std::vector<ScoringOptions::alignment_t> alignment;
	ScoringOptions::Orientation orientation;
	ScoringOptions::ScoreFunc score_func;

	void parse_alignment(const std::string& alignment_str) {
		std::istringstream ss(alignment_str);
		std::string token;
		std::set<ScoringOptions::alignment_t> unique_alignments;

		int res = 0;
		while (std::getline(ss, token, ',')) {
			int num = std::stoi(token);
			unique_alignments.insert(static_cast<ScoringOptions::alignment_t>(abs(num)));
		}

		std::copy(unique_alignments.begin(), unique_alignments.end(), std::back_inserter(alignment));
	}

	Options(int argc, char** argv) : args(argc, argv) {
		using std::string;

		auto&& positional = args.positional();
		if (positional.empty()) print_usage_and_exit();
		fasta_path = fs::path(positional[0].data());

		basename = args.get<string>("basename", (fasta_path.parent_path() / fasta_path.stem()).string());

		int max_heptad_displacement = args.get<int>("max-heptad-displacement", 0);
		auto alignment_str = args.get<string>("alignment");
		if (alignment_str) {
			if (max_heptad_displacement != 0) {
				std::cout << "alignment and max_heptad_displacement can not be specified at the same time\n";
				exit(1);
			}
			parse_alignment(alignment_str.value());
		}
		else {
			for (int i = 0; i <= max_heptad_displacement; i++) {
				alignment.push_back(+7 * i);
			}
		}

		string orientation_str = args.get<string>("orientation", "parallel");
		if (orientation_str == "parallel") {
			orientation = ScoringOptions::Orientation::parallel;
		}
		else if (orientation_str == "antiparallel") {
			orientation = ScoringOptions::Orientation::antiparallel;
		}
		else if (orientation_str == "both") {
			orientation = ScoringOptions::Orientation::both;
		}
		else {
			print_usage_and_exit();
		}

		auto score_func_str = args.get<string>("score-func", "potapov");
		if (score_func_str == "potapov") {
			score_func = ScoringOptions::ScoreFunc::potapov;
		}
		else if (score_func_str == "bcipa") {
			score_func = ScoringOptions::ScoreFunc::bcipa;
		}
		else {
			print_usage_and_exit();
		}
	}
};