#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <chrono>
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <thread>
#include <vector>
#include <omp.h>
#include <string>
#include "Interaction.h"

#include "options.h"

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef __unix
#define O_BINARY 0
#include <unistd.h>
#include <sys/mman.h>
#else
#include "mman.h"
#include <io.h>
#endif

#ifdef _MSC_VER
#define open _open
#define creat _creat
#define lseek64 _lseeki64
#define O_CREAT _O_CREAT
#define O_BINARY _O_BINARY
#define O_RDWR _O_RDWR
#define O_TRUNC _O_TRUNC
#define S_IWRITE _S_IWRITE
#endif

#define MAXN 33000

using namespace std;


void check(void);
size_t N, num_of_threads;
string fasta[MAXN];
string names[MAXN];
string fasta_name;
Interaction *inter;
int read_fasta(void);
void calc(int, int, int);
std::string get_name(int);

typedef std::chrono::system_clock clock_type;
std::chrono::time_point<clock_type> get_clock()
{
	return clock_type::now();
}

float time_diff(std::chrono::time_point<clock_type> from, std::chrono::time_point<clock_type> to)
{
	std::chrono::duration<float> dur = to - from;
	return dur.count();
}

void* mmap_wrapper(const char* fname, size_t size)
{
	int fd = open(fname, O_CREAT | O_BINARY | O_RDWR | O_TRUNC, 0644);
	if (fd == -1)
	{
		printf("ERROR open()-ing file:%s\n", strerror(errno));
		exit(1);
	}
	lseek64(fd, size - 1, SEEK_SET);
	write(fd, "", 1);
	return mmap(nullptr, size, PROT_WRITE, MAP_SHARED, fd, 0);
}

std::vector<std::string> split(std::string text, char sep) {
	std::vector<std::string> tokens;
	std::size_t start = 0, end = 0;
	while ((end = text.find(sep, start)) != std::string::npos) {
		if (end != start) {
			tokens.push_back(text.substr(start, end - start));
		}
		start = end + 1;
	}
	if (end != start) {
		tokens.push_back(text.substr(start));
	}
	return tokens;
}

int main(int argc, char **argv) {
	auto start = get_clock();

	if (argc > 1)
	{
		fasta_name = string(argv[1]);
		options::parse(argc - 1, argv + 1);
	}
	else
	{
		options::usage(argv);
		exit(0);
	}

	N = 0;
	int ret = read_fasta();
	if (ret < 0) {
		return 0;
	}
	//cout << N << endl;
	inter = new Interaction();
	inter->init_complete_score();

	string out_name, align_name, basename, orientation_name;

	if (fasta_name.length() >= 5 && 0 == fasta_name.compare(fasta_name.length() - 6, 6, ".fasta"))
	{
		basename = fasta_name.substr(0, fasta_name.length() - 6);
	}
	else
	{
		basename = fasta_name;
	}

	char orientation = options::get("orientation", 'B');
	bool check_parallel, check_antiparallel;
	switch (orientation) {
	case 'B':
		check_parallel = check_antiparallel = true; 
		break;
	case 'P':
		check_parallel = true;
		check_antiparallel = false;
		break;
	case 'A':
		check_parallel = false;
		check_antiparallel = true;
		break;
	default:
		cerr << "Unrecognized option " << orientation << " for orientation. It must be P, A, or B" << endl;
		return 0;
	}

	basename = options::get("out-name", basename);
	out_name = basename + ".bin";
	align_name = basename + ".align.bin";
	orientation_name = basename + ".orientation.bin";

	auto align_str = split(options::get("align", std::string("0")), ',');
	vector<int> alignments;
	transform(align_str.begin(), align_str.end(), back_inserter(alignments), [](string s) { return stoi(s); });
	cerr << "Will test the following alignments: ";
	for (int a : alignments) cerr << a << " ";
	cerr << endl;

	size_t sz = (size_t)N*N * sizeof(float);

	float *score = (float*)mmap_wrapper(out_name.c_str(), sz);
	int8_t *align = (int8_t*)mmap_wrapper(align_name.c_str(), N*N);

	uint8_t *antiparallel = (uint8_t*)mmap_wrapper(orientation_name.c_str(), N*N);

	if (score == MAP_FAILED)
	{
		printf("ERROR mmap()-ing file:%s\n", strerror(errno));
		exit(1);
	}

	float elapsed, speed;
	size_t done = 0, prev = 0;
	auto prev_clock = get_clock();
#pragma omp parallel for  schedule(dynamic, 10)
	for (int p = 0; p < N; p++)
	{
		if (omp_get_thread_num() == 0)
		{
			elapsed = time_diff(prev_clock, get_clock());
			speed = (done - prev) / elapsed;
			printf("%zd / %zd, %.2f per sec\n", done, N, speed);
			prev = done;
			prev_clock = get_clock();
		}

		for (int q = p; q < N; q++)
		{
			pair<float, int> score_parallel = make_pair(numeric_limits<float>::infinity(), 0);
			pair<float, int> score_antiparallel = make_pair(numeric_limits<float>::infinity(), 0); 

			float scr = numeric_limits<float>::infinity();
			int alignment = 0;
			bool antiparallel_is_best = false;

			if (check_parallel) {
				score_parallel = inter->score_alignments(fasta[p], fasta[q], alignments);
				scr = score_parallel.first;
				alignment = score_parallel.second;
			}

			if (check_antiparallel) {
				string reversed(fasta[q].length(), '\0');
				reverse_copy(fasta[q].begin(), fasta[q].end(), reversed.begin());
				score_antiparallel = inter->score_alignments(fasta[p], reversed, alignments);

				if (score_antiparallel.first < scr) {
					scr = score_antiparallel.first;
					alignment = score_antiparallel.second;
					antiparallel_is_best = true;
				}
			}

			size_t idx_pq = N*p + q;
			size_t idx_qp = N*q + p;

			score[idx_pq] = scr;
			score[idx_qp] = scr;
			align[idx_pq] = alignment;
			align[idx_qp] = alignment;
			if (check_antiparallel) {
				antiparallel[idx_pq] = antiparallel_is_best;
				antiparallel[idx_qp] = antiparallel_is_best;
			}
		}
#pragma omp atomic
		done++;
	}

	munmap(score, sz);
	munmap(align, N*N);
	munmap(antiparallel, N*N);

	auto stop = get_clock();
	printf("%.2fs elapsed\n", time_diff(start, stop));

	return 0;
}


int read_fasta() {
	N = 0;
	string tmp;
	ifstream in(fasta_name.c_str());
	if (in.is_open()) {
		while (in >> tmp) {
			if (tmp.at(0) == '>') {
				names[N] = tmp.substr(1);
			}
			else {
				fasta[N] = tmp;
				N++;
			}
		}
		in.close();
		return 0;
	}
	else {
		return -1;
	}
}

string get_name(int idx) {
	string ret = "part";
	ret.append(to_string(idx));
	return ret;
}

//popravi granice i zdruzi fajlove

void calc(int from, int to, int idx) {
	ofstream out(get_name(idx).c_str());
	int p = from / N;
	int q = from%N;
	int missed = p*(p + 1) / 2;
	q += missed;
	p += q / N;

	for (int i = from; i < to; i++) {
		if (i == (p + 1)*N - p*(p + 1) / 2) {
			p++; q = p;
		}
		out << "P" << p + 1 << "," << "P" << q + 1 << " " << inter->score_complete(fasta[p], fasta[q], 0) << endl;
		q++;
	}
	out.close();
}

void check() {
	ofstream out("complete.out");
	if (!out) return;
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			out << fasta[i] << "," << fasta[j] << inter->score_complete(fasta[i], fasta[j], 0) << endl;
		}
	}
	out.close();
}
