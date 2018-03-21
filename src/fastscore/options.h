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
    --basename=PATH                        specify output base name
    --max-heptad-displacement=NUM          try shifting the peptides left and right by up to this many heptads
    --alignment=LIST                       try these alignments
    --truncate={0, 1}                      truncate the chains when aligning them, false by default
    --orientation={parallel, antiparallel, both}
    --score-func={potapov, bcipa, qcipa	   choose scoring function
				  icipa_core_vert, icipa_nter_core}          
)");
		exit(1);
	}

	fs::path fasta_path, current_executable_path;
	std::string basename;
	std::vector<ScoringOptions::alignment_t> alignment;
	ScoringOptions::Orientation orientation;
	ScoringOptions::ScoreFunc score_func;
	bool truncate;

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

		basename = args.get<string>(
			"basename", 
			args.get<string>("out-name", 
							(fasta_path.parent_path() / fasta_path.stem()).string()));

		int max_heptad_displacement = args.get<int>("max-heptad-displacement", 0);
		auto alignment_str = args.get<string>("alignment");
		if(!alignment_str) alignment_str = args.get<string>("align");
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

		truncate = args.get<bool>("truncate", false);

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
		else if (score_func_str == "qcipa") {
			score_func = ScoringOptions::ScoreFunc::qcipa;
		}
		else if (score_func_str == "icipa_core_vert") {
			score_func = ScoringOptions::ScoreFunc::icipa_core_vert;
		}
		else if (score_func_str == "icipa_nter_core") {
			score_func = ScoringOptions::ScoreFunc::icipa_nter_core;
		}
		else {
			print_usage_and_exit();
		}
	}

	void print_parsed() {
		using std::cout;
		using std::endl;
		cout << "Fasta path is " << fasta_path.string() << endl;
		cout << "Output basename is " << basename << endl;

		cout << "Will test the following alignments:  ";
		for_each(alignment.rbegin(), alignment.rend(), [&](int a) {cout << -a << " "; });
		bool skip_zero = (alignment[0] == 0);
		for_each(alignment.begin() + (skip_zero?1:0), alignment.end(), [&](int a) {cout << a << " "; });
		cout << endl;

		cout << "Considered orientations are " << [](ScoringOptions::Orientation o) {
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

		cout << "Chains are " << (truncate ? "" : "not ") << "truncated\n";

		cout << "Used scoring function is " << [](ScoringOptions::ScoreFunc sfunc) {
			switch (sfunc) {
			case ScoringOptions::ScoreFunc::potapov:
				return "potapov";
			case ScoringOptions::ScoreFunc::bcipa:
				return "bcipa";
			case ScoringOptions::ScoreFunc::qcipa:
				return "qcipa";
			case ScoringOptions::ScoreFunc::icipa_core_vert:
				return "icipa_core_vert";
			case ScoringOptions::ScoreFunc::icipa_nter_core:
				return "icipa_nter_core";
			default:
				return "invalid";
			}
		}(score_func) << endl;
	}
};