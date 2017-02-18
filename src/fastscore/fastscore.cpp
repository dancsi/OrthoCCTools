#include <iostream>
#include <atomic>
#include <chrono>

#include "common/experimental_cxx_features.h"

#include "flags.h"

#include "ScoringHelper.h"
#include "ScoringEnginePotapov.h"
#include "ScoringEngineBCIPA.h"
#include "common/PeptideSet.h"
#include "common/InteractionMatrix.h"

using namespace std;
namespace fs = std::experimental::filesystem;

void print_usage_and_exit() {
	puts(R"(
USAGE: fastscore INPUT [OPTIONS...]
Available options:
    --basename=PATH                 specify output base name
    --max-heptad-displacement=NUM   try shifting the peptides left and right by this many heptads
    --orientation={parallel, antiparallel, both}
	--score-func={potapov, bcipa}   choose scoring function
)");
	exit(1);
}

template<typename ScoringEngineType>
void score_pairs(PeptideSet& ps, ScoringEngineType& sc, int max_heptad_displacement, ScoringOptions::Orientation orientation, InteractionMatrix& im) {
	auto n = ps.size();
	int items_to_process = n*(n + 1) / 2;
	int items_processed_total = 0, items_processed_prev = 0;
	int reporting_interval = items_to_process / 100;

	auto start = chrono::high_resolution_clock::now();
	auto time_prev = start;
	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			auto score = sc.score(ps[i].sequence, ps[j].sequence, max_heptad_displacement, orientation);
			im[i][j] = im[j][i] = score.score;
			items_processed_total++;
		}

		if (items_processed_total - items_processed_prev > reporting_interval) {
			auto time_current = chrono::high_resolution_clock::now();
			chrono::duration<double> time_elapsed = time_current - time_prev;
			float speed = 100.0 * (items_processed_total - items_processed_prev) / time_elapsed.count() / items_to_process;

			printf("%.2f %%, %.2f %% / sec\n", 100.0 * items_processed_total / items_to_process, speed);

			items_processed_prev = items_processed_total;
			time_prev = time_current;
		}
	}

	auto stop = chrono::high_resolution_clock::now();

	chrono::duration<double> duration = stop - start;
	cout << duration.count() << endl;
}

int main(int argc, char **argv) {
	const flags::args args(argc, argv);

	auto&& positional = args.positional();
	if (positional.empty()) print_usage_and_exit();
	auto fasta_path = fs::path(positional[0].data());

	string basename{ (fasta_path.parent_path() / fasta_path.stem()).string() };
	basename = args.get<std::string>("basename").value_or(basename);

	int max_heptad_displacement = args.get<int>("max-heptad-displacement", 0);

	string orientation_str = args.get<string>("orientation", "parallel");
	auto orientation = ScoringOptions::Orientation::invalid;
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

	PeptideSet ps(fasta_path.string());
	InteractionMatrix im{ basename + ".bin", ps.size() };

	auto score_func = args.get<string>("score-func", "potapov");
	if (score_func == "potapov") {
		ScoringHelper<ScoringEnginePotapov> sc;
		score_pairs(ps, sc, max_heptad_displacement, orientation, im);
	}
	else if (score_func == "bcipa") {
		ScoringHelper<ScoringEngineBCIPA> sc;
		score_pairs(ps, sc, max_heptad_displacement, orientation, im);
	}
	else {
		print_usage_and_exit();
	}

	return 0;
}