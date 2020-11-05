#include <stdlib.h>

#include "scoring/ScoringHelper.h"

#include "scoring/ScoringEngineBCIPA.h"
#include "scoring/ScoringEngineQCIPA.h"
#include "scoring/ScoringEngineICIPA.h"
#include "scoring/ScoringEnginePotapov.h"

extern "C" int __cxa_atexit(void (*func) (void *), void * arg, void * dso_handle);
extern "C" {
	int __cxa_thread_atexit(void (*func) (void *), void * arg, void * dso_handle) {
		return __cxa_atexit(func, arg, dso_handle);
	}
}

[[clang::export_name("alloc_string")]] char* alloc_string(int len) {
	char* ptr = new char[len+1];
	return ptr;
}

[[clang::export_name("free_string")]] void free_string(char* ptr) {
	delete[] ptr;
}

template<typename ScoringEngineType>
float score_func_wrapper(char *s1, int n1, char *s2, int n2, int align, int truncate, int orientation) {
	static ScoringHelper<ScoringEngineType> helper;	
	return helper.score(
		s1, s2, 
		{align},
		truncate,
		static_cast<ScoringOptions::Orientation>(orientation)
	).score;
}

#define DEFINE_FUNC(func_name, class_name) \
	[[clang::export_name(#func_name)]] float func_name(char *s1, int n1, char *s2, int n2, int align, int truncate, int orientation) { \
		return score_func_wrapper<class_name>(s1, n1, s2, n2, align, truncate, orientation); \
	}

DEFINE_FUNC(bcipa, ScoringEngineBCIPA)
DEFINE_FUNC(icipa_core_vert, ScoringEngineICIPACoreVert)
DEFINE_FUNC(icipa_nter_core, ScoringEngineICIPANterCore)