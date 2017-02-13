#pragma once

#include <map>
#include <vector>
#include <string>

#define TRACE_LEVEL 1
#define TRACE_MASK (TRACE_MASK_CLIQUE|TRACE_MASK_INITIAL)
#define PEPTIDE_LENGTH 200

extern std::vector<std::pair<int, int> > vertices, initial_set;
extern std::vector<std::string> reverse_id;

extern std::string initial_set_fname;
extern std::map<std::string, int> peptide_id;

int get_id(std::string peptide);