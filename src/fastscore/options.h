#pragma once

#include "flags.h"
#include "scoring/ScoringHelper.h"

#include <algorithm>
#include <cctype>
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

	fs::path fasta_path, current_executable_path;
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

		current_executable_path = argv[0];

		auto&& positional = args.positional();
		if (positional.empty()) print_usage_and_exit();

		fasta_path = fs::path(positional[0].data());
		if (!fs::exists(fasta_path)) {
			std::cout << "The specified fasta path "<<fasta_path.string()<<" does not exist\n";
			exit(1);
		}

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
		char first_char = std::toupper(orientation_str[0]);
		switch (first_char) {
		case 'P':
			orientation = ScoringOptions::Orientation::parallel; break;
		case 'A':
			orientation = ScoringOptions::Orientation::antiparallel; break;
		case 'B':
			orientation = ScoringOptions::Orientation::both; break;
		default:
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

	void print_parsed() {
		using std::cerr;
		using std::endl;
		cerr << "Fasta path is " << fasta_path.string() << endl;
		cerr << "Output basename is " << basename << endl;

		cerr << "Considered alignments (absolute values)\n  ";
		for (auto x : alignment) cerr << x << " ";
		cerr << endl;

		cerr << "Considered orientations are " << [](ScoringOptions::Orientation o) {
			switch (o) {
			case ScoringOptions::Orientation::parallel:
				return "parallel";
			case ScoringOptions::Orientation::antiparallel:
				return "antiparallel";
			case ScoringOptions::Orientation::both:
				return "both";
			default:
				return "invalid";
			}
		}(orientation) << endl;

		cerr << "Used scoring function is " << [](ScoringOptions::ScoreFunc sfunc) {
			switch (sfunc) {
			case ScoringOptions::ScoreFunc::potapov:
				return "potapov";
			case ScoringOptions::ScoreFunc::bcipa:
				return "bcipa";
			default:
				return "invalid";
			}
		}(score_func) << endl;
	}
};