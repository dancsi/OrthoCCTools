#pragma once

#include <format>
#include <string>

#include "common/PeptideSet.h"
#include "common/SpecialMatrices.h"
#include "scoring/ScoringHelper.h"

enum class OutputFormat
{
	binary, csv
};

using ScoringOptions::Orientation;
using ScoringOptions::alignment_t;
using OrientationMatrix = SquareMatrix<ScoringOptions::Orientation>;
using AlignmentMatrix = SquareMatrix<alignment_t>;

std::ostream& operator<<(std::ostream& os, Orientation orientation) {
	os << [orientation] {
		switch (orientation) {
		case Orientation::parallel:
			return "P";
		case Orientation::antiparallel:
			return "A";
		case Orientation::both:
			return "B";
		case Orientation::invalid:
			return "I";
		}
	}();
	return os;
}

void save_csv(
	const std::string& basename,
	const PeptideSet& peptideSet,
	const InteractionMatrix& interactionMatrix,
	const OrientationMatrix& orientationMatrix,
	const AlignmentMatrix& alignmentMatrix) 
{
	size_t n = peptideSet.size();

	std::string outputName = basename + ".csv";
	std::ofstream outputFile(outputName);

	outputFile << "ID1,ID2,score,orientation,alignment\n";
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j <= i; j++) {
			outputFile << peptideSet[i].id << ',' << peptideSet[j].id << ',' << interactionMatrix[i][j] << ',' << orientationMatrix[i][j] << ',' << alignmentMatrix[i][j] << '\n';
		}
	}
}

void save_bin(
	const std::string& basename,
	const PeptideSet& peptideSet,
	const InteractionMatrix& interactionMatrix,
	const OrientationMatrix& orientationMatrix,
	const AlignmentMatrix& alignmentMatrix) 
{
	interactionMatrix.to_binary(basename + ".bin");
	orientationMatrix.to_binary(basename + ".orientation.bin");
	alignmentMatrix.to_binary(basename + ".align.bin");
}

void save(OutputFormat outputFormat,
	const std::string& basename,
	const PeptideSet& peptideSet,
	const InteractionMatrix& interactionMatrix,
	const OrientationMatrix& orientationMatrix,
	const AlignmentMatrix& alignmentMatrix) {
	switch (outputFormat) {
		case OutputFormat::binary:
			save_bin(basename, peptideSet, interactionMatrix, orientationMatrix, alignmentMatrix);
			break;
		case OutputFormat::csv:
			save_csv(basename, peptideSet, interactionMatrix, orientationMatrix, alignmentMatrix);
			break;
	}
}