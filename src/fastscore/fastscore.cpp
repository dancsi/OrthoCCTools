#include <iostream>
#include <atomic>
#include <chrono>

#include "common/experimental_cxx_features.h"

#include "flags.h"

#include "ScoringEngine.h"
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
)");
	exit(1);
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
	auto orientation = ScoringEngine::Orientation::invalid;
	if (orientation_str == "parallel") {
		orientation = ScoringEngine::Orientation::parallel;
	}
	else if (orientation_str == "antiparallel") {
		orientation = ScoringEngine::Orientation::antiparallel;
	}
	else if (orientation_str == "both") {
		orientation = ScoringEngine::Orientation::both;
	}
	else {
		print_usage_and_exit();
	}

	ScoringEngine sc;

	PeptideSet ps(fasta_path.string());
	auto n = ps.size();
	
	InteractionMatrix im{ basename + ".bin", ps.size() };

	auto start = chrono::high_resolution_clock::now();
	//#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			auto score = sc.score(ps[i].sequence, ps[j].sequence, max_heptad_displacement, orientation);
			im[i][j] = score.score;
		}
	}
	auto stop = chrono::high_resolution_clock::now();
	chrono::duration<double> duration = stop - start;
	cout << duration.count() << endl;
	
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			cout << ps[i].id << "," << ps[j].id << "," << im[i][j] << endl;
		}
	}

	return 0;
}