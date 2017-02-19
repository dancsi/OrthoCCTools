#pragma once

#include "flags.h"
#include "ScoringHelper.h"

#include <string>

namespace fs = std::experimental::filesystem;

struct Options {
	const flags::args args;

	void print_usage_and_exit() {
		puts(R"(
USAGE: fastscore INPUT [OPTIONS...]
Available options:
    --basename=PATH                 specify output base name
    --max-heptad-displacement=NUM   try shifting the peptides left and right by this many heptads
	--alignment=-7*N,-7*(N-1),...,-7,0,7,...,7*(N-1),7*N
    --orientation={parallel, antiparallel, both}
	--score-func={potapov, bcipa}   choose scoring function
)");
		exit(1);
	}

	fs::path fasta_path;
	std::string basename;
	int max_heptad_displacement;
	ScoringOptions::Orientation orientation;
	ScoringOptions::ScoreFunc score_func;

	int parse_alignment(const std::string& alignment) {
		std::istringstream ss(alignment);
		std::string token;

		int res = 0;
		while (std::getline(ss, token, ',')) {
			int num = std::stoi(token);
			res = std::max(res, abs(num) / 7);
		}

		return res;
	}

	Options(int argc, char** argv) : args(argc, argv) {
		using std::string;

		auto&& positional = args.positional();
		if (positional.empty()) print_usage_and_exit();
		fasta_path = fs::path(positional[0].data());

		basename = args.get<string>("basename", (fasta_path.parent_path() / fasta_path.stem()).string());

		max_heptad_displacement = args.get<int>("max-heptad-displacement", 0);
		auto alignment = args.get<string>("alignment");
		if (alignment) {
			max_heptad_displacement = parse_alignment(alignment.value());
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