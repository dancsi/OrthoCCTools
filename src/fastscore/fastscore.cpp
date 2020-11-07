#include <iostream>
#include <atomic>
#include <chrono>

#include "options.h"

#include "scoring/ScoringHelper.h"
#include "scoring/ScoringEnginePotapov.h"
#include "scoring/ScoringEngineBCIPA.h"
#include "scoring/ScoringEngineQCIPA.h"
#include "scoring/ScoringEngineICIPA.h"

#include "common/PeptideSet.h"
#include "common/SpecialMatrices.h"

using namespace std;
namespace fs = std::filesystem;

using ScoringOptions::alignment_t;
using OrientationMatrix = SquareMatrix<ScoringOptions::Orientation>;
using AlignmentMatrix = SquareMatrix<alignment_t>;

template<typename ScoringEngineType>
void score_pairs(
	PeptideSet& ps,
	ScoringHelper<ScoringEngineType>& sc,
	std::vector<alignment_t>& alignment,
	bool truncate,
	ScoringOptions::Orientation orientation,
	InteractionMatrix& im,
	OrientationMatrix& om,
	AlignmentMatrix& am
) {
	auto n = ps.size();
	int items_to_process = n * (n + 1) / 2;
	int items_processed_total = 0, items_processed_prev = 0;
	int reporting_interval = items_to_process / 100;

	auto start = chrono::high_resolution_clock::now();
	auto time_prev = start;
#pragma omp parallel for schedule(dynamic, 4)
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			auto score = sc.score(ps[i].sequence, ps[j].sequence, alignment, truncate, orientation);
			im[i][j] = im[j][i] = score.score;
			om[i][j] = om[j][i] = score.orientation;
			am[i][j] = score.alignment; am[j][i] = -score.alignment;
			items_processed_total++;
		}
#pragma omp critical 
		{
			if (items_processed_total - items_processed_prev > reporting_interval) {
				auto time_current = chrono::high_resolution_clock::now();
				chrono::duration<double> time_elapsed = time_current - time_prev;
				float speed = 100.0 * (items_processed_total - items_processed_prev) / time_elapsed.count() / items_to_process;

				fprintf(stderr, "%.2f %%, %.2f %% / sec\n", 100.0 * items_processed_total / items_to_process, speed);

				items_processed_prev = items_processed_total;
				time_prev = time_current;
			}
		}
	}

	auto stop = chrono::high_resolution_clock::now();

	chrono::duration<double> duration = stop - start;
	cout << duration.count() << endl;
}

int main(int argc, char **argv) {
	Options options(argc, argv);
	options.print_parsed();

	auto& alignment = options.alignment;
	auto truncate = options.truncate;
	auto orientation = options.orientation;

	auto fasta_path = options.fasta_path.string();
	auto basename = options.basename;

	PeptideSet ps(fasta_path);
	InteractionMatrix im{ basename + ".bin", ps.size() };
	OrientationMatrix om{ basename + ".orientation.bin", ps.size() };
	AlignmentMatrix am{ basename + ".align.bin", ps.size() };

	auto score_func = options.score_func;
	if (score_func == ScoringOptions::ScoreFunc::potapov) {
		ScoringHelper<ScoringEnginePotapov> sc(
			ScoringEnginePotapov::find_weights_file(
				fs::path(argv[0]).parent_path()
			)
		);
		score_pairs(ps, sc, alignment, truncate, orientation, im, om, am);
	}
	else if (score_func == ScoringOptions::ScoreFunc::bcipa) {
		ScoringHelper<ScoringEngineBCIPA> sc;
		score_pairs(ps, sc, alignment, truncate, orientation, im, om, am);
	}
	else if (score_func == ScoringOptions::ScoreFunc::qcipa) {
		ScoringHelper<ScoringEngineQCIPA> sc;
		score_pairs(ps, sc, alignment, truncate, orientation, im, om, am);
	}
	else if (score_func == ScoringOptions::ScoreFunc::icipa_core_vert) {
		ScoringHelper<ScoringEngineICIPACoreVert> sc;
		score_pairs(ps, sc, alignment, truncate, orientation, im, om, am);
	}
	else if (score_func == ScoringOptions::ScoreFunc::icipa_nter_core) {
		ScoringHelper<ScoringEngineICIPANterCore> sc;
		score_pairs(ps, sc, alignment, truncate, orientation, im, om, am);
	}

	return 0;
}