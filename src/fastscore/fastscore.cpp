#include <algorithm>
#include <execution>
#include <iostream>
#include <atomic>
#include <chrono>
#include <thread>

#include "io.h"
#include "options.h"

#include "scoring/ScoringHelper.h"
#include "scoring/ScoringEnginePotapov.h"
#include "scoring/ScoringEngineBCIPA.h"
#include "scoring/ScoringEngineQCIPA.h"
#include "scoring/ScoringEngineICIPA.h"

#include "common/ParallelFor.h"
#include "common/PeptideSet.h"
#include "common/SpecialMatrices.h"
#include "common/TriangluarIndex.h"

using namespace std;

template<typename ScoringEngineType>
void score_pairs(
	PeptideSet& ps,
	vector<alignment_t>& alignment,
	bool truncate,
	ScoringOptions::Orientation orientation,
	InteractionMatrix& im,
	OrientationMatrix& om,
	AlignmentMatrix& am
) {
	ScoringHelper<ScoringEngineType> sc{};
	auto n = ps.size();
	int items_to_process = n * (n + 1) / 2;
	int items_processed_total = 0, items_processed_prev = 0;
	int reporting_interval = items_to_process / 100;

	auto start = chrono::high_resolution_clock::now();
	auto time_prev = start;

	auto indices = TriangularIndices(n);
	parallel_for(indices.begin(), indices.end(), [&](auto idx) 
		{
			auto [i, j] = idx;
			auto score = sc.score(ps[i].sequence, ps[j].sequence, alignment, truncate, orientation);
			im[i][j] = im[j][i] = score.score;
			om[i][j] = om[j][i] = score.orientation;
			am[i][j] = score.alignment; am[j][i] = -score.alignment;
		}
	);

	auto stop = chrono::high_resolution_clock::now();

	cout << "Done in " << chrono::duration_cast<chrono::seconds>(stop - start).count() << " seconds" << endl;
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
	InteractionMatrix im(ps.size());
	OrientationMatrix om(ps.size());
	AlignmentMatrix am(ps.size());

	switch (options.score_func)
	{
	case ScoringOptions::ScoreFunc::potapov:
		score_pairs<ScoringEnginePotapov>(ps, alignment, truncate, orientation, im, om, am);
		break;
	case ScoringOptions::ScoreFunc::bcipa:
		score_pairs<ScoringEngineBCIPA>(ps, alignment, truncate, orientation, im, om, am);
		break;
	case ScoringOptions::ScoreFunc::qcipa:
		score_pairs<ScoringEngineQCIPA>(ps, alignment, truncate, orientation, im, om, am);
		break;
	case ScoringOptions::ScoreFunc::icipa_core_vert:
		score_pairs<ScoringEngineICIPACoreVert>(ps, alignment, truncate, orientation, im, om, am);
		break;
	case ScoringOptions::ScoreFunc::icipa_nter_core:
		score_pairs<ScoringEngineICIPANterCore>(ps, alignment, truncate, orientation, im, om, am);
		break;
	}

	save(options.output_format, basename, ps, im, om, am);

	return 0;
}