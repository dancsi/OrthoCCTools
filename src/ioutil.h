#pragma once

#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "mcqd_para/MaximumCliqueBase.h"
#include "mcqd_para/ParallelMaximumClique.h"
#include "mcqd_para/BB_GreedyColorSort.h"
#include "mcqd_para/McrBB.h"

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
#define fstat64 _fstat64
#define O_CREAT _O_CREAT
#define O_BINARY _O_BINARY
#define O_RDWR _O_RDWR
#define O_TRUNC _O_TRUNC
#define S_IWRITE _S_IWRITE
#endif

bool parse_input(FILE* fin, char* id1, char* id2, float& score)
{
	char buf[2 * PEPTIDE_LENGTH + 10];
	char *ret = fgets(buf, sizeof(buf) - 1, fin);
	if (feof(fin) || !ret) return false;

	int len = strlen(buf);
	if (len<5) return false;

	int i = 0;
	while (i <len &&buf[i] != ',')
	{
		*id1++ = buf[i];
		i++;
	}
	if (i == len) return false;
	*id1 = 0; i++;
	while (i<len && buf[i] != ',')
	{
		*id2++ = buf[i];
		i++;
	}
	if (i == len) return false;
	*id2 = 0; i++;
	score = atof(buf + i);
	return true;
}

void dump_dimacs(std::vector<std::vector<char>>& conn, const char *fname)
{
	FILE* fout = fopen(fname, "w");
	int n = conn.size(), m = conn[0].size();
	fprintf(fout, "p edge %d %d\n", n, m);
	for (int i = 0;i<n;i++)
	{
		for (int j = 0;j<i;j++)
		{
			if (conn[i][j])
				fprintf(fout, "e %d %d\n", i + 1, j + 1);
		}
	}
	fclose(fout);
}

void print_clique(const string out_name, const BitstringSet& clique, Graph<BitstringSet>& graph)
{
	stringstream ss;
	bool comma = false;
	BitstringSet bs = clique;

	graph.remap(bs); //VERY IMPORTANT!!!

	while (bs.size()>0)
	{
		int i = bs.nextSetBit();
		bs.remove(i);
		{
			auto p = vertices[i];
			ss << reverse_id[p.first] << "," << reverse_id[p.second] << "\n";
		}
	}

	for (pair<int, int>& p : initial_set)
	{
		ss << reverse_id[p.first] << "," << reverse_id[p.second] << "\n";
	}

	printf("%s", ss.str().c_str());

	FILE* fout = fopen(out_name.c_str(), "w");
	fprintf(fout, "%s", ss.str().c_str());
	fclose(fout);
}

void read_initial_set()
{
	if (initial_set_fname.empty()) return;
	FILE* fin = fopen(initial_set_fname.c_str(), "r");
	if (fin == nullptr) return;
	fprintf(stderr, "Reading initial set from %s\n", initial_set_fname.c_str());

	char p1[PEPTIDE_LENGTH], p2[PEPTIDE_LENGTH];
	while (fscanf(fin, "%[^,],%s\n", p1, p2) != EOF)
	{
		auto it1 = peptide_id.find(string(p1));
		auto it2 = peptide_id.find(string(p2));
		if (it1 == peptide_id.end() || it2 == peptide_id.end())
		{
			initial_set.clear();
			fprintf(stderr, "Invalid peptide ID %s, aborting.", (it1 == peptide_id.end()) ? p1 : p2);
			exit(1);
		}

		auto newpair = make_pair(it1->second, it2->second);
		if (newpair.first > newpair.second) swap(newpair.first, newpair.second);
		initial_set.push_back(newpair);
	}

	fclose(fin);
}

float** allocate_matrix_pointers(float* data, int n)
{
	float** mat = new float*[n];
	for (int i = 0;i < n;i++) mat[i] = data + i*n;
	return mat;
}

float **read_scores_binary(std::string score_file)
{
	int fd = open(score_file.c_str(), O_RDONLY, 0600);

	if (fd == -1)
	{
		perror("Error opening file for writing");
		exit(EXIT_FAILURE);
	}

	struct _stat64 fileInfo = { 0 };

	if (fstat64(fd, &fileInfo) == -1)
	{
		perror("Error getting the file size");
		exit(EXIT_FAILURE);
	}

	if (fileInfo.st_size == 0)
	{
		fprintf(stderr, "Error: File is empty, nothing to do\n");
		exit(EXIT_FAILURE);
	}

	printf("File size is %ji\n", (intmax_t)fileInfo.st_size);

	float* score_storage = (float*)mmap(0, fileInfo.st_size, PROT_READ, MAP_SHARED, fd, 0);
	if (score_storage == MAP_FAILED)
	{
		close(fd);
		perror("Error mmapping the file");
		exit(EXIT_FAILURE);
	}
	size_t n_peptides = (size_t)sqrt(fileInfo.st_size / 4);
	char buf[20];
	for (int i = 0;i < n_peptides;i++)
	{
		sprintf(buf, "P%d", i+1);
		get_id(buf);
	}
	return allocate_matrix_pointers(score_storage, n_peptides);
}

float **read_scores_plaintext(std::string score_file)
{
	float *score_storage = new float[4100 * 4100];
	float **score = allocate_matrix_pointers(score_storage, 4100);

	FILE* fin = fopen(score_file.c_str(), "r");
	int id1, id2;
	float _score;
	char id1s[PEPTIDE_LENGTH], id2s[PEPTIDE_LENGTH];

	while (parse_input(fin, id1s, id2s, _score))
	{
		id1 = get_id(id1s);
		id2 = get_id(id2s);

		score[id1][id2] = score[id2][id1] = _score;
	}
	fclose(fin);
	return score;
}

float **read_scores(std::string score_file)
{
	if (score_file.find(".bin") != std::string::npos)
	{
		return read_scores_binary(score_file);
	}
	else
	{
		return read_scores_plaintext(score_file);
	}
}