#include <iostream>
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

#ifdef _MSC_VER
#include "mman.h"
#include <io.h>
#define open _open
#define creat _creat
#define O_CREAT _O_CREAT
#define O_BINARY _O_BINARY
#define O_RDWR _O_RDWR
#define O_TRUNC _O_TRUNC
#define S_IWRITE _S_IWRITE
#else
#include <sys/mman.h>
#endif

#ifdef __unix
#define O_BINARY 0
#include <unistd.h>
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

int main(int argc, char **argv) {
	clock_t start = clock();
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

	float scr, elapsed, speed;
	size_t done = 0, prev = 0;
	clock_t prev_clock = clock();
	const float alpha = 0.5;
#pragma omp parallel for private(scr)
	for (int p = 0;p < N;p++)
	{
		if(omp_get_thread_num() == 0)
		{
			elapsed = (float)(clock() - prev_clock) / CLOCKS_PER_SEC;
			speed = (done - prev) / elapsed;
			printf("%zd / %zd, %.2f per sec\n", done, N, speed);
			prev = done;
			prev_clock = clock();
		}

		for (int q = p;q < N;q++)
		{
			scr = inter->score_complete(fasta[p], fasta[q]);
			score[N*p + q] = scr;
			score[N*q + p] = scr;
		}
#pragma omp atomic
			done++;
	}

	munmap(score, sz);
	clock_t stop = clock();
	printf("%.2fs elapsed\n", float(stop-start) / CLOCKS_PER_SEC);

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
