#include <iostream>
#include <chrono>
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <thread>
#include <omp.h>
#include "Interaction.h"

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

typedef std::chrono::high_resolution_clock clock_type;
std::chrono::time_point<clock_type> get_clock()
{
	return clock_type::now();
}

float time_diff(std::chrono::time_point<clock_type> from, std::chrono::time_point<clock_type> to)
{
	std::chrono::duration<float> dur = to - from;
	return dur.count();
}

int main(int argc, char **argv) {
	auto start = get_clock();
	if (argc > 1) {
		fasta_name = string(argv[1]);
	}
	N = 0;
	int ret = read_fasta();
	if (ret < 0) {
		return 0;
	}
	//cout << N << endl;
	inter = new Interaction("alll.in");
	inter->init_complete_score();

	int fd = open("output.bin", O_CREAT | O_BINARY | O_RDWR | O_TRUNC, 0644);
	if(fd == -1)
	{
		printf("ERROR open()-ing file:%s\n", strerror(errno));
		exit(1);
	}

	size_t sz = (size_t)N*N*sizeof(float);
	lseek64(fd, sz, SEEK_SET);
	write(fd, "", 1);

	float *score = (float*)mmap(nullptr, sz, PROT_WRITE, MAP_SHARED, fd, 0);
	if(score == MAP_FAILED)
	{
		printf("ERROR mmap()-ing file:%s\n", strerror(errno));
		exit(1);
	}

	float elapsed, speed;
	size_t done = 0, prev = 0;
	auto prev_clock = get_clock();
	const float alpha = 0.5;
#pragma omp parallel for  schedule(guided)
	for (int p = 0;p < N;p++)
	{
		if(omp_get_thread_num() == 0)
		{
			elapsed = time_diff(prev_clock, get_clock());
			speed = (done - prev) / elapsed;
			printf("%zd / %zd, %.2f per sec\n", done, N, speed);
			prev = done;
			prev_clock = get_clock();
		}

		for (int q = p;q < N;q++)
		{
			float scr = inter->score_complete(fasta[p], fasta[q]);
			score[N*p + q] = scr;
			score[N*q + p] = scr;
		}
#pragma omp atomic
			done++;
	}

	munmap(score, sz);
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

	for (int i = from;i < to;i++) {
		if (i == (p + 1)*N - p*(p + 1) / 2) {
			p++; q = p;
		}
		out << "P" << p + 1 << "," << "P" << q + 1 << " " << inter->score_complete(fasta[p], fasta[q]) << endl;
		q++;
	}
	out.close();
}

void check() {
	ofstream out("complete.out");
	if (!out) return;
	for (int i = 0;i < N;i++) {
		for (int j = i;j < N;j++) {
			out << fasta[i] << "," << fasta[j] << inter->score_complete(fasta[i], fasta[j]) << endl;
		}
	}
	out.close();
}
