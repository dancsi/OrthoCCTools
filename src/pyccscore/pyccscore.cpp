#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "scoring/ScoringHelper.h"

#include "scoring/ScoringEngineBCIPA.h"
#include "scoring/ScoringEngineQCIPA.h"
#include "scoring/ScoringEngineICIPA.h"
#include "scoring/ScoringEnginePotapov.h"

namespace py = pybind11;
namespace fs = std::filesystem;

template<typename ScoringEngineType>
void register_scoring_engine(pybind11::module& module, std::shared_ptr<ScoringHelper<ScoringEngineType>> helper, std::string_view name) {
	using namespace py::literals;
	using namespace ScoringOptions;

	module.def(name.data(), [helper](
		std::string_view chain1,
		std::string_view chain2,
		std::vector<alignment_t> alignment,
		bool truncate,
		Orientation orientation) {
			auto ret = helper->score(chain1, chain2, alignment, truncate, orientation);
			return std::make_tuple(ret.score, ret.alignment, ret.orientation);
		},
			"chain1"_a, "chain2"_a,
			py::arg_v("alignment", std::vector<alignment_t>{0}, "[0]"),
			"truncate"_a = false,
			py::arg_v("orientation", Orientation::parallel));
}

PYBIND11_MODULE(pyccscore, m)
{
	m.doc() = "Python bindings for coiled-coil scoring functions";

	py::enum_<ScoringOptions::Orientation>(m, "Orientation", "Peptide orientation")
		.value("PARALLEL", ScoringOptions::Orientation::parallel)
		.value("ANTIPARALLEL", ScoringOptions::Orientation::antiparallel)
		.value("BOTH", ScoringOptions::Orientation::both)
		.value("INVALID", ScoringOptions::Orientation::invalid)
		.export_values();

	auto potapov = std::make_shared<ScoringHelper<ScoringEnginePotapov>>();
	auto bcipa = std::make_shared<ScoringHelper<ScoringEngineBCIPA>>();
	auto qcipa = std::make_shared<ScoringHelper<ScoringEngineQCIPA>>();
	auto icipa_core_vert = std::make_shared<ScoringHelper<ScoringEngineICIPACoreVert>>();
	auto icipa_nter_core = std::make_shared<ScoringHelper<ScoringEngineICIPANterCore>>();

	register_scoring_engine(m, potapov, "score_potapov");
	register_scoring_engine(m, bcipa, "score_bcipa");
	register_scoring_engine(m, qcipa, "score_qcipa");
	register_scoring_engine(m, icipa_core_vert, "icipa_core_vert");
	register_scoring_engine(m, icipa_nter_core, "icipa_nter_core");
}