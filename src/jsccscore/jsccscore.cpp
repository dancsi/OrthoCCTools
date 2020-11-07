#include <stdlib.h>
#include <sstream>

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

extern "C" void fill_score_dat_buffer(char* buf); 
std::string score_dat_buffer() {
	char *buf = new char[300000];
	fill_score_dat_buffer(buf);
	std::string str{buf};
	delete[] buf;
	return str;
}

char aligned_oriented_score_buf[sizeof(ScoringOptions::aligned_oriented_score_t)];

template<typename T>
char* score_func_wrapper(ScoringHelper<T>& helper, char *s1, int n1, char *s2, int n2, int max_align, int truncate, int orientation) {
	std::vector<ScoringOptions::alignment_t> align;
	for(int i=0;i<=max_align;i+=7) align.push_back(i); 

	auto score = helper.score(
		s1, s2, 
		align,
		truncate,
		static_cast<ScoringOptions::Orientation>(orientation)
	);
	memcpy(aligned_oriented_score_buf, &score, sizeof(score));

	return aligned_oriented_score_buf;
}

#define DEFINE_FUNC(func_name, class_name, engine_init) \
	[[clang::export_name(#func_name)]] char* func_name(char *s1, int n1, char *s2, int n2, int max_align, int truncate, int orientation) { \
		static auto helper = ScoringHelper<class_name>(engine_init); \
		return score_func_wrapper(helper, s1, n1, s2, n2, max_align, truncate, orientation); \
	}

DEFINE_FUNC(bcipa, ScoringEngineBCIPA,)
DEFINE_FUNC(icipa_core_vert, ScoringEngineICIPACoreVert,)
DEFINE_FUNC(icipa_nter_core, ScoringEngineICIPANterCore,)
DEFINE_FUNC(potapov, ScoringEnginePotapov, score_dat_buffer())